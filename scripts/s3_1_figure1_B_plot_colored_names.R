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

# lunge_merged = readRDS(here("results/s2_1_lunge_merged.RDS"))
lunge_merged = readRDS(here("../results/s507_4_lunge_merged.RDS"))

# RENAME CLUSTERS ----



clusterorder = c('Alveolar Macrophages', 'Proliferating Alveolar Macrophages', 'Macrophages/ Monocytes', 'Mast Cells', 'Dendritic Cells', 'T Cells', 'NK Cells', 'Proliferating NK/T', 'B Cells', 'AT1', 'AT2', 'Ciliated Cells', 'Club Cells', 'Endothelial', 'ly Endothelial', 'Fibroblasts', 'Perivalscular' )


clusternames_new_pre = c("C 0" = "T Cells", 'C 1' = "Endothelial", 'C 2' = "Alveolar Macrophages", 'C 3' = "NK Cells", 'C 4' = "Macrophages/ Monocytes", 'C 5' = "Endothelial", 'C 6' = "Fibroblasts", 'C 7' = "AT2", 'C 8' = "Perivascular", 'C 9' = "B Cells", 'C 10' = "Macrophages/ Monocytes", 'C 11' = "Mast Cells", 'C 12' = "Macrophages/ Monocytes", 'C 13' = "Endothelial", 'C 14' = "AT1", 'C 15' = "Alveolar Macrophages", 'C 16' = "Ciliated Cells", 'C 17' = "Endothelial", 'C 18' = "Club Cells", "C 19" = "ly Endothelial", "C 20" = "Dendritic Cells", "C 21" = "Endothelial", "C 22" = "T Cells", "C 23" = "Proliferating NK/T", "C 24" = "Proliferating Alveolar Macrophages", "C 25" = "Endothelial", "C 26" = "T Cells", "C 27" = "B Cells")

clusternames_new = data.table(old_id  = names(clusternames_new_pre),
                              new_id = clusternames_new_pre)

venn2(clusternames_new$old_id, paste0("C ", lunge_merged$SCT_snn_res.0.4))

lunge_merged$cluster_seurat_v2 = clusternames_new[match_hk(paste0("C ", lunge_merged$SCT_snn_res.0.4), clusternames_new$old_id), new_id]

Idents(lunge_merged) = lunge_merged$cluster_seurat_v2

lunge_merged$spezies = ifelse(lunge_merged$species=="hamster", "Hamster",
                              ifelse(lunge_merged$species=="mouse", "Mouse",
                                     ifelse(lunge_merged$species=="human no FACS", "HumanK",
                                            ifelse(lunge_merged$species=="human FACS", "HumanT",lunge_merged$species))))


lunge_merged$spezies = factor(lunge_merged$spezies, levels = c("Mouse", "Hamster", "HumanK", "HumanT"))

unique(data.table(lunge_merged$spezies , lunge_merged$species))


Nogpalette2 <-  c("#CB769E", "#F7C548", "#F97E44", "#62C370", "#FB3640", "#B7245C", "#4F6D7A", "#0D3B66", "#B2675E", "#3E2F5B", "#A85C85", "#7C6A0A", "#246A73", "#644536", "#DE639A", "#5CC1BC", "#368F8B")

farben_df = data.table(tissue = sort(unique(clusternames_new$new_id)), Nogpalette2 = Nogpalette2)

clusternames_new$farben =farben_df[match_hk(clusternames_new$new_id, farben_df$tissue), Nogpalette2]

    showNA(lunge_merged@meta.data, showAllNoNA = F)
p5b =
  DimPlot(lunge_merged, reduction = "umap",pt.size = 0.8, split.by ="spezies", group.by = "cluster_seurat_v2", label = T, label.size = 3,repel = T,  cols = c(Nogpalette2)) + coord_fixed(ratio=1) +NoLegend()
p5b


require(ggrepel)
head(p5b$data)

p5b_data_pre = data.table(p5b$data, keep.rownames = T)


p5b_data_pre[, median_x := median(UMAP_1), .(spezies, cluster_seurat_v2)]
p5b_data_pre[, median_y := median(UMAP_2), .(spezies, cluster_seurat_v2)]
p5b_data_label = unique(p5b_data_pre[,.(median_x,median_y, cluster_seurat_v2, spezies)])

p5c = ggplot(p5b_data_pre, aes(UMAP_1, UMAP_2, col = cluster_seurat_v2)) + geom_point(size = 1, alpha = 0.3) + facet_grid(.~spezies) +
  geom_text_repel(data = p5b_data_label, aes(x = median_x, y = median_y, label = str_wrap(cluster_seurat_v2, width = 20)), seed = 1902, size = 4.5) +
  theme_cowplot() +
  guides(color = "none") +
  scale_color_manual(values = Nogpalette2) +
  geom_text_repel(data = p5b_data_label, aes(x = median_x, y = median_y, label = str_wrap(cluster_seurat_v2, width = 20)), seed = 1902, alpha = 0.4, col = "black", size = 4.5) +
  Seurat:::FacetTheme() +
  coord_fixed(ratio=1) +
  xlab("UMAP 1") +
  ylab("UMAP 2")

p5c



# savePlotInMultiFormats(corefilename = here("results/s508_5_ERJ_methods_v2_colored"),
#                        weite = 18,
#                        hoehe = 5,
#                        plotecode = "print(p5c)"
# )
#
#
# plot(1:nrow(clusternames_new), 1:nrow(clusternames_new), pch = " ",
#         col = clusternames_new$farben)
#
# text(1:nrow(clusternames_new), 1:nrow(clusternames_new),
#        clusternames_new$new_id, col = clusternames_new$farben)
#
# write.delim(clusternames_new, here("results/s508_5_zuordnung_ERJ_method_clustV2.txt"))
#
# readr::write_rds(lunge_merged, file = here("results/s508_5_lunge_merged_clusternames_v2.RDS"))


finalizeSkript()
