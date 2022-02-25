
table(Errint$dataset)
Errint = Errint[,Errint$dataset != 'humanSS2']
table(Errint$dataset)

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

dev.off()



# # FINALIZE SCRIPT----
finalizeSkript()
