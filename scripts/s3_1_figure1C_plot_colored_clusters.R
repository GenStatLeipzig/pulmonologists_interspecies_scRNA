# # LOAD DATA ----
#
Errint <- readRDS(".../Data/s603_2_seurat_lunge_merged_21-12-11.RDS") # processed integrated data, available at https://www.health-atlas.de/studies/54


# # CREATE DATA ----

NogpaletteApha <-  c("#CB769E", "#F7C548", "#F97E44", "#62C370", "#FB3640", "#B7245C", "#4F6D7A", "#0D3B66", "#B2675E", "#3E2F5B", "#A85C85", "#7C6A0A", "#0081AF", "#246A73", "#644536", "#DE639A", "#5CC1BC", "#368F8B")


p5b =
  DimPlot(Errint, reduction = "umap",pt.size = 0.8, split.by ="species", group.by = "celltype", label = T, label.size = 3,repel = T,  cols = c(NogpaletteApha)) + coord_fixed(ratio=1) +NoLegend()
p5b

# # COLOR LABELS----
p5b_data_pre = data.table(p5b$data, keep.rownames = T)


p5b_data_pre[, median_x := median(UMAP_1), .(species, celltype)]
p5b_data_pre[, median_y := median(UMAP_2), .(species, celltype)]
p5b_data_label = unique(p5b_data_pre[,.(median_x,median_y, celltype, species)])

p5c = ggplot(p5b_data_pre, aes(UMAP_1, UMAP_2, col = celltype)) + geom_point(size = 1, alpha = 0.3) + facet_grid(.~species) +
  geom_text_repel(data = p5b_data_label, aes(x = median_x, y = median_y, label = str_wrap(celltype, width = 20)), seed = 1902, size = 4.5) +
  theme_cowplot() +
  guides(color = "none") +
  scale_color_manual(values = NogpaletteApha) +
  geom_text_repel(data = p5b_data_label, aes(x = median_x, y = median_y, label = str_wrap(celltype, width = 20)), seed = 1902, alpha = 0.2, col = "black", size = 4.5) +
  Seurat:::FacetTheme() +
  coord_fixed(ratio=1) +
  xlab("UMAP 1") +
  ylab("UMAP 2")

p5c

dev.off()



# FINALIZE SCRIPT----
finalizeSkript()
