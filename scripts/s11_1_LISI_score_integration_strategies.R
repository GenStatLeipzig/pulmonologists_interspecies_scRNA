# # Initialize Skript ----
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
require(lisi)
require(dplyr)

# # LOAD and annotate seurat objects----
seurat_fn = here("results/s603_2_seurat_lunge_merged_21-12-11.RDS")
stopifnot(file.exists(seurat_fn))
seurat_fn



Errint_harmony =  readRDS(seurat_fn)
Errint_harmony@active.ident = factor(Errint_harmony$SCT_snn_res.0.3)

Errint_harmony$celltype <- paste("C", Errint_harmony$SCT_snn_res.0.3)
Errint_harmony@active.ident = factor(Errint_harmony$celltype)

Errint_harmony <- RenameIdents (Errint_harmony, "C 0" = "T Cells", 'C 1' = "Endothelial", 'C 2' = "Alveolar Macrophages", 'C 3' = "Macrophages/ Monocytes", 'C 4' = "NK Cells", 'C 5' = "AT2", 'C 6' = "Endothelial", 'C 7' = "Fibroblasts", 'C 8' = "B Cells", 'C 9' = "Perivascular", 'C 10' = "AT1", 'C 11' = "Macrophages/ Monocytes", 'C 12' = "Mast Cells", 'C 13' = "Ciliated Cells", 'C 14' = "Club Cells", 'C 15' = "Dendritic Cells", 'C 16' = "Macrophages/ Monocytes", 'C 17' = "ly Endothelial", 'C 18' = "Proliferating NK/T", "C 19" = "Proliferating AM", "C 20" = "Neutrophils", "C 21" = "Endothelial", "C 22" = "Fibroblasts")



reihenfolge = c('Alveolar Macrophages', 'Proliferating AM', 'Macrophages/ Monocytes', "Neutrophils", 'Dendritic Cells', 'Mast Cells', 'T Cells', 'NK Cells', 'Proliferating NK/T', 'B Cells', 'AT1', 'AT2', 'Ciliated Cells', 'Club Cells', 'Endothelial', 'ly Endothelial', 'Fibroblasts', 'Perivascular')

NogpaletteReihe <-  c("#CB769E", "#DE639A", "#A85C85", "#0081AF", "#4F6D7A", "#7C6A0A", "#368F8B", "#246A73", "#5CC1BC", "#62C370", "#F7C548", "#F97E44", "#FB3640", "#B7245C", "#0D3B66", "#3E2F5B", "#B2675E", "#644536")

venn2(Errint_harmony@active.ident, reihenfolge)
Errint_harmony@active.ident <- factor(x = Errint_harmony@active.ident, levels = reihenfolge)




p_harmony = UMAPPlot(Errint_harmony, label = T, cols = c(NogpaletteReihe), label.size = 3.0) + NoLegend()  + ggtitle("Integration via Harmony")+ xlab("UMAP 1") + ylab("UMAP 2") + scale_color_manual(drop=FALSE, values = NogpaletteReihe)
p_harmony


## alternatively processed object 
Errint_seurat = readRDS(here("results/s630_1_seurat_lunge_merged_noMITOadjSEURATAnchors_22-01-02.RDS"))

Errint_seurat_anno = Errint_seurat@meta.data %>% as.data.table()


## plot UMAP----
venn2(Errint_seurat$rev1_celltype_completed, reihenfolge)
Errint_seurat@active.ident = factor(Errint_seurat$rev1_celltype_completed, levels  = reihenfolge)

p_seurat = UMAPPlot(Errint_seurat, label = T, cols = c(NogpaletteReihe), label.size = 3.0) + NoLegend() + ggtitle("Integration via Seurat integrateData()") + xlab("UMAP 1") + ylab("UMAP 2") + scale_color_manual(drop=FALSE, values = NogpaletteReihe)
p_seurat

require(patchwork)
p_harmony + p_seurat

# LISI harmony ----
# measure for integration
stopifnot(identical(rownames(p_harmony$data), rownames(Errint_harmony@meta.data)))
stopifnot(identical(rownames(p_harmony$data), names(Errint_harmony@active.ident)))

metadata_harmony = cbind(p_harmony$data[, c("UMAP_1","UMAP_2")], celltype = Errint_harmony@active.ident , run10x = Errint_harmony@meta.data$run10x)
head(metadata_harmony)
uniqueN(metadata_harmony$run10x, na.rm = T)
uniqueN(metadata_harmony$celltype, na.rm = T)

lisi_harmony <- compute_lisi(X = metadata_harmony[,c("UMAP_1","UMAP_2")], 
                             meta_data = metadata_harmony, 
                             label_colnames = c("celltype","run10x")
)
                              

sapply(lisi_harmony, median)


# LISI seurat ----
stopifnot(identical(rownames(p_seurat$data), rownames(Errint_seurat@meta.data)))
stopifnot(identical(rownames(p_seurat$data), names(Errint_seurat@active.ident)))

metadata_seurat = cbind(p_seurat$data[, c("UMAP_1","UMAP_2")], celltype = Errint_seurat@active.ident , run10x = Errint_seurat@meta.data$run10x)
head(metadata_seurat)
uniqueN(metadata_seurat$run10x, na.rm = T)
uniqueN(metadata_seurat$celltype, na.rm = T)

lisi_seurat <- compute_lisi(X = metadata_seurat[,c("UMAP_1","UMAP_2")], 
                             meta_data = metadata_seurat, 
                             label_colnames = c("celltype","run10x")
)


sapply(lisi_seurat, median)





## plot run10x-----
Errint_harmony@active.ident = factor(Errint_harmony$run10x)

p_harmony2 = UMAPPlot(Errint_harmony, label = F) + NoLegend() + ggtitle("Integration via harmony") + xlab("UMAP 1") + ylab("UMAP 2") 
p_harmony2


Errint_seurat@active.ident = factor(Errint_seurat$run10x)

p_seurat2 = UMAPPlot(Errint_seurat, label = F) + NoLegend() + ggtitle("Integration via Seurat integrateData()") + xlab("UMAP 1") + ylab("UMAP 2") 
p_seurat2


## plot LISI dets-----
Errint_harmony$lisi_run10x = lisi_harmony$run10x 

p_harmony3 = FeaturePlot(Errint_harmony, features = "lisi_run10x", cols = c("blue", 'orange'))  + ggtitle("Integration via harmony") + xlab("UMAP 1") + ylab("UMAP 2") 
p_harmony3


Errint_seurat$lisi_run10x = lisi_seurat$run10x 

p_seurat3 = FeaturePlot(Errint_seurat, features = "lisi_run10x", cols = c("blue", 'orange'))  + ggtitle("Integration via seurat") + xlab("UMAP 1") + ylab("UMAP 2") 
p_seurat3



## joint plot -----
### density ----
plotdat = rbind(lisi_harmony %>% mutate(method = "harmony", celltypename = metadata_harmony$celltype),
                lisi_seurat %>% mutate(method = "seurat" , celltypename = metadata_seurat$celltype)
                ) %>% data.table()
ggplot(plotdat, aes(run10x, col  = method)) + geom_density() + theme_minimal()


### comparison median, global and celltype level ----

plotdat1 = plotdat[, .(iLISI = median(run10x)), .(method)]
plotdat1


p3 = ggplot(plotdat1, aes(reorder(method, -iLISI), iLISI, fill = method)) + geom_col(position = "dodge") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) + scale_y_continuous(breaks = pretty_breaks(10))

plotdat2 = plotdat[, .(iLISI = median(run10x)), .(celltypename, method)]
plotdat2



p4 = ggplot(plotdat2, aes(reorder(celltypename, -iLISI), iLISI, fill = method)) + geom_col(position = "dodge") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) + scale_y_continuous(breaks = pretty_breaks(10))

p3+p4 + patchwork::plot_layout(widths = c(1,4),guides = "collect" ) & xlab("")





finalizeSkript()



