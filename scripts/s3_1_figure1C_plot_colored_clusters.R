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


# # create PLOTDATA----
Idents(lunge_merged) = lunge_merged$cluster_seurat_v2


Nogpalette3 <-  c("#CB769E", "#DE639A", "#A85C85", "#7C6A0A", "#4F6D7A", "#368F8B",
                  "#246A73", "#5CC1BC", "#62C370", "#F7C548", "#F97E44", "#FB3640",
                  "#B7245C", "#0D3B66", "#3E2F5B", "#B2675E", "#644536")


p5b =
  DimPlot(lunge_merged, reduction = "umap",pt.size = 0.8, split.by ="species", group.by = "cluster_seurat_v2", label = T, label.size = 3,repel = T,  cols = c(Nogpalette3)) + coord_fixed(ratio=1) +NoLegend()
p5b

# # COLOR LABELS----
p5b_data_pre = data.table(p5b$data, keep.rownames = T)


p5b_data_pre[, median_x := median(UMAP_1), .(species, cluster_seurat_v2)]
p5b_data_pre[, median_y := median(UMAP_2), .(species, cluster_seurat_v2)]
p5b_data_label = unique(p5b_data_pre[,.(median_x,median_y, cluster_seurat_v2, species)])

p5c = ggplot(p5b_data_pre, aes(UMAP_1, UMAP_2, col = cluster_seurat_v2)) + geom_point(size = 1, alpha = 0.3) + facet_grid(.~species) +
  geom_text_repel(data = p5b_data_label, aes(x = median_x, y = median_y, label = str_wrap(cluster_seurat_v2, width = 20)), seed = 1902, size = 4.5) +
  theme_cowplot() +
  guides(color = "none") +
  scale_color_manual(values = Nogpalette3) +
  geom_text_repel(data = p5b_data_label, aes(x = median_x, y = median_y, label = str_wrap(cluster_seurat_v2, width = 20)), seed = 1902, alpha = 0.2, col = "black", size = 4.5) +
  Seurat:::FacetTheme() +
  coord_fixed(ratio=1) +
  xlab("UMAP 1") +
  ylab("UMAP 2")

p5c


jpeg(here("results/s3_1_integrated_cluters.jpeg"),width =  18,height =  5,unit = 'in', quality = 100, res = 150)
print(p5c)
dev.off()



# FINALIZE SCRIPT----
finalizeSkript()
