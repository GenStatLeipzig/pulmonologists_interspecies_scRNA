rm(list = ls())

library(SingleCellSignalR)
library(Seurat)
library(generics)
library(here)
require(toolboxH)
require(stringr)
require(ggplot2)
require(ggthemes)
require(scales)
require(cowplot)

# # LOAD seurat objekt laden----
seurat_fn = here("results/s603_2_seurat_lunge_merged_21-12-11.RDS")
stopifnot(file.exists(seurat_fn))
seurat_fn



Errint_harmony =  readRDS(seurat_fn)
Errint_harmony@active.ident = factor(Errint_harmony$SCT_snn_res.0.3)

Errint_harmony$celltype <- paste("C", Errint_harmony$SCT_snn_res.0.3)
Errint_harmony@active.ident = factor(Errint_harmony$celltype)

Errint_harmony <- RenameIdents (Errint_harmony, "C 0" = "T Cells", 'C 1' = "Endothelial", 'C 2' = "Alveolar Macrophages", 'C 3' = "Macrophages/ Monocytes", 'C 4' = "NK Cells", 'C 5' = "AT2", 'C 6' = "Endothelial", 'C 7' = "Fibroblasts", 'C 8' = "B Cells", 'C 9' = "Perivascular", 'C 10' = "AT1", 'C 11' = "Macrophages/ Monocytes", 'C 12' = "Mast Cells", 'C 13' = "Ciliated Cells", 'C 14' = "Club Cells", 'C 15' = "Dendritic Cells", 'C 16' = "Macrophages/ Monocytes", 'C 17' = "ly Endothelial", 'C 18' = "Proliferating NK/T", "C 19" = "Proliferating AM", "C 20" = "Neutrophils", "C 21" = "Endothelial", "C 22" = "Fibroblasts")



Errint_harmony$celltype <- factor(Errint_harmony@active.ident, levels = sort(unique(Errint_harmony@active.ident)))

nogpalette3 = c("#CB769E", "#DE639A", "#A85C85", "#0081AF", "#4F6D7A", "#7C6A0A", "#368F8B", "#246A73", "#5CC1BC", "#62C370", "#F7C548", "#F97E44", "#FB3640", "#B7245C", "#0D3B66", "#3E2F5B", "#B2675E", "#644536")


p_harmony = DimPlot(Errint_harmony, label = T)  + guides(color = "none") + scale_color_manual(values =nogpalette3 ) + ggtitle("Integration via R-package Harmony")+ xlab("UMAP 1") + ylab("UMAP 2")
p_harmony


# alternative laden
Errint_seurat = readRDS(here("results/s630_1_seurat_lunge_merged_noMITOadjSEURATAnchors_22-01-02.RDS"))
Errint_harmony$celltype <- factor(Errint_seurat$rev1_celltype_completed, levels = sort(unique(Errint_seurat$rev1_celltype_completed)))

p_seurat = DimPlot(Errint_seurat, label = T,group.by = "rev1_celltype_completed", repel = T) + guides(color = "none") + scale_color_manual(values =nogpalette3 ) +   scale_fill_manual(values = nogpalette3) + ggtitle("Integration via R-package Seurat integrate(Data)") + xlab("UMAP 1") + ylab("UMAP 2")
p_seurat

require(patchwork)
p_harmony + p_seurat


# jpeg(here("results/s608d_1_supplFig2_dotplotMarker.jpeg"), 14,18, units = "in", quality = 100, res = 300)
# p2
# dev.off()
# 
toolboxH::savePlotInMultiFormats(here("results/s631_1_supplFig_vgl_UMAP_harmony_seuratAnchor"), weite = 16,8, plotecode = 'plot(p_harmony + p_seurat)')



finalizeSkript()
