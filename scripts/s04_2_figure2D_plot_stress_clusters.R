
rm(list = ls())
require(toolboxH)
require(here)
require(Seurat)
require(ggplot2)
require(ggthemes)
require(dplyr)
require(tidyr)
require(cowplot)
require(paletteer)
require(patchwork)


Errint = readRDS(here("results/s603_2_seurat_lunge_merged_21-12-11.RDS"))

table(Errint$dataset)
Errint = Errint[,Errint$dataset != 'humanSS2']
table(Errint$dataset)

Errint@active.ident = factor(Errint$SCT_snn_res.0.3)

Errint$celltype <- paste("C", Errint$SCT_snn_res.0.3)
Errint@active.ident = factor(Errint$celltype)

Errint <- RenameIdents (Errint, "C 0" = "T Cells", 'C 1' = "Endothelial", 'C 2' = "Alveolar Macrophages", 'C 3' = "Macrophages/ Monocytes", 'C 4' = "NK Cells", 'C 5' = "AT2", 'C 6' = "Endothelial", 'C 7' = "Fibroblasts", 'C 8' = "B Cells", 'C 9' = "Perivascular", 'C 10' = "AT1", 'C 11' = "Macrophages/ Monocytes", 'C 12' = "Mast Cells", 'C 13' = "Ciliated Cells", 'C 14' = "Club Cells", 'C 15' = "Dendritic Cells", 'C 16' = "Macrophages/ Monocytes", 'C 17' = "ly Endothelial", 'C 18' = "Proliferating NK/T", "C 19" = "Proliferating AM", "C 20" = "Neutrophils", "C 21" = "Endothelial", "C 22" = "Fibroblasts")

Errint$celltype <- paste(Errint@active.ident)
table(Errint$SCT_snn_res.0.3, Errint@active.ident)

unique(Errint@meta.data[,c("SCT_snn_res.0.3", "celltype")])
reihenfolge = c('Alveolar Macrophages', 'Proliferating AM', 'Macrophages/ Monocytes', "Neutrophils", 'Dendritic Cells', 'Mast Cells', 'T Cells', 'NK Cells', 'Proliferating NK/T', 'B Cells', 'AT1', 'AT2', 'Ciliated Cells', 'Club Cells', 'Endothelial', 'ly Endothelial', 'Fibroblasts', 'Perivascular')
NogpaletteReihe <-  c("#CB769E", "#DE639A", "#A85C85", "#0081AF", "#4F6D7A", "#7C6A0A", "#368F8B", "#246A73", "#5CC1BC", "#62C370", "#F7C548", "#F97E44", "#FB3640", "#B7245C", "#0D3B66", "#3E2F5B", "#B2675E", "#644536")
my_levels <- c(reihenfolge)
Errint@active.ident <- factor(x = Errint@active.ident, levels = my_levels)

stressgene <- c("HSPA8", "FOS", "DUSP1", "IER3", "EGR1", "FOSB", "HSPB1", "ATF3")

p1 = FeaturePlot(Errint, features = stressgene, pt.size = 0.2,
                 ncol = 1, split.by = "species", keep.scale = "feature" ,cols = c("blue", "yellow"))
p1

Idents(Errint)  = Errint$species
stressgene_Errint_species = AverageExpression(Errint, features = stressgene)
stressgene_Errint_species$SCT

Errint <- AddModuleScore(
  object = Errint,
  features = list(stressgene),
  name = 'Stressgenes',
  assay = "SCT"
) # https://github.com/satijalab/seurat/issues/3022
# wird komischerweise nicht als Stressgenes, sondern als Stressgenes1 angefuegt

Errint.meta.data = data.table(Errint@meta.data, keep.rownames = T)
Errint.meta.data[,n_spezies := uniqueN(species), celltype]
ggplot(Errint.meta.data, aes(species, Stressgenes1, fill = species)) + geom_violin()


p2 = DotPlot(Errint, features = c("Stressgenes1", stressgene),group.by = "celltype" , split.by = "species", cols = c("blue", 'blue', 'blue', "blue","blue", 'blue'))  + RotatedAxis()
p2


p2$data$species = str_split(p2$data$id, "_") %>% sapply(., "[", 2)
p2$data$tissue = str_split(p2$data$id, "_") %>% sapply(., "[", 1)

newdat = data.table(p2$data)
newdat
newdat[,species := factor(species, levels = c("hamster", "human", "monkey",  "mouse",  "pig", "rat"))]
newdat[,genes := ifelse(features.plot =="Stressgenes1", "stress-\nscore", as.character(features.plot))]
newdat[,genes := factor(genes,levels = c("stress-\nscore", setdiff(sort(unique(genes)), "stress-\nscore")))]


reihenfolge = c('Alveolar Macrophages', 'Proliferating AM', 'Macrophages/ Monocytes', "Neutrophils", 'Dendritic Cells', 'Mast Cells', 'T Cells', 'NK Cells', 'Proliferating NK/T', 'B Cells', 'AT1', 'AT2', 'Ciliated Cells', 'Club Cells', 'Endothelial', 'ly Endothelial', 'Fibroblasts', 'Perivascular')

venn2(reihenfolge, newdat$tissue)

newdat[,tissue :=factor(tissue, levels = reihenfolge)]
p3 = ggplot(newdat, aes(x = tissue, y = species)) + geom_point(aes(size = pct.exp, col = (avg.exp.scaled)))+ coord_flip() + facet_grid(.~genes) + theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.4)) + scale_color_paletteer_c("viridis::viridis") + labs(col = "scaled\nexpression", size = "%\nexpressed") + ylab("") + xlab("")
p3

ggsave(here("results/s608b_1_Stress_dotplots2_10xonly.pdf"), width = 16, height = 6)
dev.off()
#
# ggsave("C:/Users/Admin/Documents/Doktorarbeit/Projekt neu/Paper Progress/ERJ Methods/ERR Update/Results/Stress_dotplots2.pdf", width = 16, height = 6)
#
# dev.off()
#


sessionInfo()
