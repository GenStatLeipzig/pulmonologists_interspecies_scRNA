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




# # CREATE STRESS SCORE ----

stressgene = c("ATF3", "BTG1", "BTG2", "DUSP1", "EGR1", "FOS", "FOSB", "HSP90AB1",
               "HSPA1A", "HSPA1B", "HSPA8", "HSPB1", "IER2", "IER3", "JUN",
               "JUNB", "JUND")


Idents(lunge_merged)  = lunge_merged$species

stressgene_lunge_merged_species = AverageExpression(lunge_merged, features = stressgene)
stressgene_lunge_merged_species$SCT

uninformativ_stressgene_pre = apply(stressgene_lunge_merged_species$SCT, 1, function(x) any(x==0))
uninformativ_stressgene = uninformativ_stressgene_pre[uninformativ_stressgene_pre==T]
colplotLikeExcel(stressgene_lunge_merged_species$SCT, myround = 3)

uninformativ_stressgene = uninformativ_stressgene_pre[uninformativ_stressgene_pre==T] %>% names
uninformativ_stressgene

informativ_stressgene = uninformativ_stressgene_pre[uninformativ_stressgene_pre==F] %>% names
informativ_stressgene

lunge_merged <- AddModuleScore(
  object = lunge_merged,
  features = list(informativ_stressgene),
  name = 'Stressgenes',
  assay = "SCT"
)

# # PLOT STRESS SCORE----

lunge_merged.meta.data = data.table(lunge_merged@meta.data, keep.rownames = T)




p5b =  DimPlot(lunge_merged, reduction = "umap",pt.size = 0.8, split.by ="species", group.by = "cluster_seurat_v2", label = T, label.size = 3,repel = T) + coord_fixed(ratio=1) +NoLegend()
p5b # grundplot clusters


head(p5b$data)

p5b_data_pre = data.table(p5b$data, keep.rownames = T)


p5b_data_pre[, median_x := median(UMAP_1), .(species, cluster_seurat_v2)]
p5b_data_pre[, median_y := median(UMAP_2), .(species, cluster_seurat_v2)]
p5b_data_label = unique(p5b_data_pre[,.(median_x,median_y, cluster_seurat_v2, species)])

p5b_data_pre[, Stressgenes1 := lunge_merged.meta.data[match_hk(p5b_data_pre$rn, lunge_merged.meta.data$rn),Stressgenes1] ]

setorder(p5b_data_pre, Stressgenes1)

p8e = ggplot(p5b_data_pre, aes(UMAP_1, UMAP_2, col =(Stressgenes1)))      + geom_point(size = 1, alpha = 1) + facet_grid(.~species)  + theme_cowplot() +
  geom_text_repel(data = p5b_data_label, aes(x = median_x, y = median_y, label = str_wrap(cluster_seurat_v2, width = 20)), seed = 1902, alpha = 0.5, col = "black") +
  Seurat:::FacetTheme() +  scale_color_gradient2(low = "blue", mid = "blue", high = "yellow", midpoint =0.5) + theme(legend.position = "top", legend.key.width = unit(1, units = "cm"))+ labs(col = "Cell-specific Stress-score") # (3+scale(Stressgenes1))^1.5)

p8e


jpeg(here("results/s4_1_stress_score.jpeg"),width =  18,height =  5,unit = 'in', quality = 100, res = 150)
print(p8e)
dev.off()



# # FINALIZE SCRIPT----
finalizeSkript()
