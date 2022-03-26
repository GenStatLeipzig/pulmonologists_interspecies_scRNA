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



Errint =  readRDS(seurat_fn)
Errint@active.ident = factor(Errint$SCT_snn_res.0.3)

Errint$celltype <- paste("C", Errint$SCT_snn_res.0.3)
Errint@active.ident = factor(Errint$celltype)

Errint <- RenameIdents (Errint, "C 0" = "T Cells", 'C 1' = "Endothelial", 'C 2' = "Alveolar Macrophages", 'C 3' = "Macrophages/ Monocytes", 'C 4' = "NK Cells", 'C 5' = "AT2", 'C 6' = "Endothelial", 'C 7' = "Fibroblasts", 'C 8' = "B Cells", 'C 9' = "Perivascular", 'C 10' = "AT1", 'C 11' = "Macrophages/ Monocytes", 'C 12' = "Mast Cells", 'C 13' = "Ciliated Cells", 'C 14' = "Club Cells", 'C 15' = "Dendritic Cells", 'C 16' = "Macrophages/ Monocytes", 'C 17' = "ly Endothelial", 'C 18' = "Proliferating NK/T", "C 19" = "Proliferating AM", "C 20" = "Neutrophils", "C 21" = "Endothelial", "C 22" = "Fibroblasts")

Errint$celltype <- paste(Errint@active.ident)


canonical_markers <- c("ADGRE1", "CD68", "MRC1", "MARCO", "CD79B", "MS4A1", "CD3E", "CD3D", "NCR1", "NKG7", "S100A8", "FLT3", "EPCAM", "RTKN2", "SFTPD", "SFTPC", "FOXJ1", "CCDC17", "PECAM1", "TMEM100", "MMRN1", "FBLN1","PDGFRA",  "ACTA2", "CNN1", "COX4I2", "CSPG4", "SCGB1A1", "MS4A2", "TOP2A")

p1 = DotPlot(Errint, features = canonical_markers, split.by = "species",group.by = "celltype" ,cols = c("red", "blue", "#8A9747", "#F37252", "#4E1339" , "black"))


p2_data = data.table(p1$data, keep.rownames = F)

p2_data[,species:=ifelse(grepl('human$',id), "human",
                         ifelse(grepl('mouse$',id), "mouse",
                                ifelse(grepl('monkey$',id), "monkey",
                                       ifelse(grepl('rat$',id), "rat",
                                              ifelse(grepl('pig$',id), "pig",
                                                     ifelse(grepl('hamster$',id), "hamster",id))))))]

p2_data[,.N,.(id, species)] %>% data.frame()

p2_data[,celltype:= str_replace(id, paste0("_",species), "")]
setnames(p2_data, "features.plot", "Markergenes")

p2 = ggplot(p2_data, aes(Markergenes, id, size =pct.exp, alpha = avg.exp.scaled, col = species))      +
  geom_point() +
  facet_grid(celltype~., scales = "free")  +
  theme_cowplot() +
  Seurat:::FacetTheme() +  scale_color_manual(values = c("red", "blue", "#8A9747", "#F37252", "#4E1339" , "black")) +
  theme(legend.position = "right",
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle =90, vjust = 0.4, hjust = 1)) +
  labs(col = "Species",
       alpha = "Average\nexpression",
       size = "Percent\nexpressed") 

p2


jpeg(here("results/s608d_1_supplFig2_dotplotMarker.jpeg"), 14,18, units = "in", quality = 100, res = 300)
p2
dev.off()

toolboxH::savePlotInMultiFormats(here("results/s608d_1_supplFig2_dotplotMarker"), 14,18, plotecode = 'plot(p2)')



finalizeSkript()
