#' ---
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'

# # INITIATE SCRIPT ----
rm(list = ls())


r_on_server = F

if(r_on_server) {
  rootpath = "/net/ifs1/san_projekte/projekte/"
  .libPaths("/net/ifs1/san_projekte/projekte/genstat/07_programme/rpackages/forostar")
  toolboxH::initializeSkript(computer = "forostar" )
} else {
  rootpath =  "R:"
  toolboxH::initializeSkript(computer = "local" )
}
require(toolboxH) # https://github.com/holgerman/toolboxH
knitr::opts_chunk$set(cache = F, results = "markup", echo=T, fig.width = 10, fig.height = 7 )

packages2load = c("magrittr","knitr", "ggplot2", "scales", "ggthemes", "assertr","Seurat", 'paletteer',"pheatmap", "cowplot", "here", "plotly", "ggrepel")

new.packages <- packages2load[!(packages2load %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

for(i in packages2load) {

  suppressPackageStartupMessages(library(i, character.only = TRUE))

}

# # LOAD DATA ----
#
lunge_merged = readRDS(here("results/s1_3_lunge_merged.RDS")) # processed integrated data, available at https://www.health-atlas.de/studies/54

lunge_merged@active.ident = factor(lunge_merged$cluster_seurat_v2)



# # PLOT CANONICAL MARKER GENES----

canonical_markers <- c("ADGRE1", "MAFB", "MARCO", "CD79A", "CD79B", "CD3E", "NKG7", "FLT3", "EPCAM", "RTKN2", "SFTPD", "FOXJ1", "PECAM1", "MMRN1", "INMT", "DCN", "ACTA2", "GUCY1A3", "SCGB1A1", "MS4A2", "TOP2A")

p1 = DotPlot(lunge_merged, features = canonical_markers, split.by = "species",group.by = "cluster_seurat_v2" ,cols = c("red", "blue", "#8A9747", "black"))


p2_data = data.table(p1$data, keep.rownames = F)

p2_data[,species:=ifelse(grepl('human charite',id), "human charite",
                         ifelse(grepl('human travaglini',id), "human travaglini",
                                ifelse(grepl('mouse',id), "mouse",
                                       ifelse(grepl('hamster',id), "hamster",id))))]

p2_data[,celltype:= str_replace(id, paste0("_",species), "")]
setnames(p2_data, "features.plot", "Markergenes")


p2 = ggplot(p2_data, aes(Markergenes, id, size =pct.exp, alpha = avg.exp.scaled, col = species))      +
  geom_point() +
  facet_grid(celltype~., scales = "free")  +
  theme_cowplot() +
  Seurat:::FacetTheme() +  scale_color_manual(values = c("red", "black",  "#8A9747","blue")) +
  theme(legend.position = "right",
        strip.text.y = element_blank(),
      axis.text.x = element_text(angle =90, vjust = 0.4, hjust = 1)) +
 labs(col = "Species",
       alpha = "Average\nexpression",
       size = "Percent\nexpressed") # (3+scale(Stressgenes1))^1.5)

p2


jpeg(here("results/s5_1_canonical_markers.jpeg"),width =  10,height =  13,unit = 'in', quality = 100, res = 150)
print(p2)
dev.off()



# # FINALIZE SCRIPT----
finalizeSkript()
