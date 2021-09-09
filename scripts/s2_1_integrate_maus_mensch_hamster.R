# # INITIIEREN ####
rm(list = ls())
#
# .libPaths("/net/ifs1/san_projekte/projekte/genstat/07_programme/rpackages/forostar3")


filename  = "s507_4_integrate_maus_mensch_hamster"

r_on_server = T


if(r_on_server) {
  # setwd( "/net/ifs1/san_projekte/projekte/genstat/02_projekte/2005_hamsterMiceHuman_scRNA_methods")
  rootpath = "/net/ifs1/san_projekte/projekte/"
  .libPaths("/net/ifs1/san_projekte/projekte/genstat/07_programme/rpackages/forostar3")
  # setwd('/net/ifs1/san_projekte/projekte/genstat/02_projekte/2005_hamsterMiceHuman_scRNA_methods/')


  toolboxH::initializeSkript(computer = "forostar3" )
} else {
  rootpath =  "R:"
  toolboxH::initializeSkript(computer = "local" )
}

require(toolboxH)

knitr::opts_chunk$set(cache = F, results = "markup", echo=T, fig.width = 10, fig.height = 7 )
#### weitere benoetigte packages und funktionen laden
# packages2load = c("AnnotationDbi","magrittr","knitr", "ggplot2", "scales", "ggthemes", "assertr","Seurat", "scater", "Matrix", "plotly", "stringr", 'SingleCellExperiment', 'scran', "scales",   'UpSetR', "scran", "batchelor",   'tidytext', 'paletteer','TabulaMurisData',   'scmap', 'gprofiler2', 'SingleR' ,'PCAtools', 'org.Mm.eg.db'  , "AUCell", 'DelayedMatrixStats', "GSEABase", "cluster", "pheatmap", "cowplot", "here")

packages2load = c("magrittr","knitr", "ggplot2", "scales", "ggthemes", "assertr","Seurat", 'paletteer',"pheatmap", "cowplot", "here", "plotly", "harmony")

for(i in packages2load) {
  # if(i %in% rownames(installed.packages()) == FALSE) {BiocManager::install(i, update = F )}
  suppressPackageStartupMessages(library(i, character.only = TRUE))
}
require(ggrepel)
require(ggalluvial)
.libPaths()

# BiocManager::install("ggalluvial")
# library(devtools)
# install_github("immunogenomics/harmony") # worked on aman






files.sources = dir(here("../../07_programme/github/scRNATexMex/R/"), full.names = T)
sapply(files.sources, source)

# PIPELINE adopted from https://github.com/immunogenomics/harmony/issues/41   kvshams commented on 11 Jun 2020 •



#laden ----
annogene = fread(here("results/s506_3_gene_annotation_orthologues_harmonized4integration.txt"))

lung_list = readr::read_rds(here("results/s506_3_list_mensch_maus_hamster_harmonized4integration.RDS"))

names(lung_list)


for (i in 1:length(lung_list)) {
  message("SCT on ", i, " i.e. ", names(lung_list)[i], "\n-------------------------------------------------------")

  lung_list[[i]] <- SCTransform(lung_list[[i]], verbose = T)
}


# sm = SCTransform(lung_list[[4]], verbose = T)

print(table(names(warnings() ))) # "iteration limit reached" ChristophH commented on 22 May 2019



## < table of extent 0 >
# These warnings are showing that there are some genes for which it is hard to reliably estimate theta (presumably because of very few non-zero observations). Usually we don't worry about these warnings too much, since we regularize the parameters in a later step, thus averaging out uncertainty of individual gene parameters. https://github.com/ChristophH/sctransform/issues/25


list.features <- SelectIntegrationFeatures(object.list = lung_list, nfeatures = 3000)
str(list.features)
##  chr [1:3000] "SFTPC" "SFTPB" "SFTPA2" "CCL2" "C1QB" "MGP" "C1QA" "NKG7" ...
lunge_merged <- merge(lung_list[[1]],
                      y = lung_list[2:length(lung_list)],
                      project = "lung",
                      merge.data = TRUE)
Assays(lunge_merged)

DefaultAssay(lunge_merged)
## [1] "SCT"
DefaultAssay(lunge_merged) = "SCT"
DefaultAssay(lunge_merged)
## [1] "SCT"
if(r_on_server==F) rm(lung_list);gc()
##              used    (Mb) gc trigger    (Mb)   max used    (Mb)
## Ncells    3206512   171.3    4910965   262.3    4910965   262.3
## Vcells 1660345703 12667.5 2810511476 21442.6 2389498827 18230.5
## PCA----
VariableFeatures(lunge_merged) <- list.features
str(VariableFeatures(lunge_merged) )
##  chr [1:3000] "SFTPC" "SFTPB" "SFTPA2" "CCL2" "C1QB" "MGP" "C1QA" "NKG7" ...
lunge_merged <- RunPCA(object = lunge_merged, assay = "SCT", npcs = 50 ,features = list.features)


table(lunge_merged@meta.data$orig.ident)

table(lunge_merged@meta.data$species)

# Todo in schritt vorher----

table(lunge_merged@meta.data$tissue, useNA = "always")

# # harmony-----
table(lunge_merged@meta.data$run10x)

lunge_merged <- RunHarmony(object = lunge_merged,
                           assay.use = "SCT",
                           reduction = "pca",
                           dims.use = 1:50,
                           group.by.vars = "run10x",
                           plot_convergence = TRUE) # test , kmeans_init_nstart=20, kmeans_init_iter_max=100 according to   https://github.com/immunogenomics/harmony/issues/25  see also https://stackoverflow.com/questions/21382681/kmeans-quick-transfer-stage-steps-exceeded-maximum



# UPDATE mit s507_4 nur noch warning Warnung: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.SCT.harmony; see ?make.names for more details on syntax validity
#
lunge_merged <- RunUMAP(object = lunge_merged, assay = "SCT", reduction = "harmony", dims = 1:50)

lunge_merged <- FindNeighbors(object = lunge_merged, assay = "SCT", reduction = "harmony", dims = 1:50)
#
lunge_merged <- FindClusters(object = lunge_merged, resolution = 0.4)

lunge_merged@active.ident %>% head

table(lunge_merged$SCT_snn_res.0.4)

# lunge_merged@meta.data$species = ifelse(grepl("^lung", lunge_merged@meta.data$orig.ident), "human no FACS",
#                                        ifelse(grepl("Human", lunge_merged@meta.data$orig.ident), "human FACS",
#                                               ifelse(grepl("ma_d0", lunge_merged@meta.data$orig.ident),"hamster",
#                                                      ifelse(grepl("X_wt", lunge_merged@meta.data$orig.ident), "mouse", NA))))
# lunge_merged$species %>% str
lunge_merged$species = factor(lunge_merged$species, levels = c("hamster", "mouse", "human no FACS", "human FACS"))
lunge_merged$species_tissue = paste0(lunge_merged$species,"_", lunge_merged$tissue)
table(lunge_merged$species, useNA = "ifany")

plotSankey(lunge_merged@meta.data[, c("SCT_snn_res.0.4", "tissue")], spalte4color = 'SCT_snn_res.0.4', gap.width = 0.5)



plotSankey(lunge_merged@meta.data[lunge_merged@meta.data$SCT_snn_res.0.4 %in% 0:10, c("SCT_snn_res.0.4", "tissue")], spalte4color = 'SCT_snn_res.0.4', gap.width = 0.5)


plotSankey(lunge_merged@meta.data[lunge_merged@meta.data$SCT_snn_res.0.4 %in% 11:20, c("SCT_snn_res.0.4", "tissue")], spalte4color = 'SCT_snn_res.0.4', gap.width = 0.5)


plotSankey(lunge_merged@meta.data[lunge_merged@meta.data$SCT_snn_res.0.4 %in% 21:30, c("SCT_snn_res.0.4", "tissue")], spalte4color = 'SCT_snn_res.0.4', gap.width = 0.5)


pdf(here("results/s507_4_cluster_sankey_zuordnung_v2.pdf"), width = 9, height = 14)
plotSankey(lunge_merged@meta.data[lunge_merged@meta.data$SCT_snn_res.0.4 %in% 0:6, c("species","SCT_snn_res.0.4", "species_tissue", "species")], spalte4color = 'SCT_snn_res.0.4', gap.width = 0.5)
plotSankey(lunge_merged@meta.data[lunge_merged@meta.data$SCT_snn_res.0.4 %in% 7:12, c("species","SCT_snn_res.0.4", "species_tissue", "species")], spalte4color = 'SCT_snn_res.0.4', gap.width = 0.5)
plotSankey(lunge_merged@meta.data[lunge_merged@meta.data$SCT_snn_res.0.4 %in% 13:18, c("species","SCT_snn_res.0.4", "species_tissue", "species")], spalte4color = 'SCT_snn_res.0.4', gap.width = 0.5)
plotSankey(lunge_merged@meta.data[lunge_merged@meta.data$SCT_snn_res.0.4 %in% 19:24, c("species","SCT_snn_res.0.4", "species_tissue", "species")], spalte4color = 'SCT_snn_res.0.4', gap.width = 0.5)
plotSankey(lunge_merged@meta.data[lunge_merged@meta.data$SCT_snn_res.0.4 %in% 25:30, c("species","SCT_snn_res.0.4", "species_tissue", "species")], spalte4color = 'SCT_snn_res.0.4', gap.width = 0.5)


dev.off()

## nochmal sankey als einer -----

sankeywhole = melt.data.table(data.table(lunge_merged@meta.data, keep.rownames = T), id.vars = c("rn" ), measure.vars =  c("species","SCT_snn_res.0.4", "species_tissue"),  variable.name = "category")


sankeywhole$Freq = 1
sankeywhole

# sankeywhole[,category := factor(category, levels = c("channel", "magnetic.selection", "free_annotation", "location")) ]

p_sankeywhole = ggplot(sankeywhole,#[rn %in% unique(rn)[1:10000]],
                       aes(x = category      , stratum = value , alluvium = rn, y = Freq, fill = value, label = value)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  # geom_text(stat = "stratum", size = 3) +
  geom_text_repel( stat = "stratum",segment.colour = "black",size = 3,min.segment.length=0) +
  theme(legend.position = "none") +
  ggtitle("Integrated data set")


p_sankeywhole



pdf(here("results/s507_4_cluster_sankey_zuordnung_v2_ueberblick.pdf"), width = 24, height = 21)
print(p_sankeywhole)

dev.off()

# autoguess clusters-----

table(lunge_merged$SCT_snn_res.0.4)

autoassign = data.table(cluster  = lunge_merged$SCT_snn_res.0.4,
                        tissue = lunge_merged$tissue)[,.N, .(cluster, tissue)]

setorder(autoassign, -N)

autoassign_v0 = autoassign[duplicated(cluster)==F][order(cluster)]
autoassign_v0

autoassign_v1 = autoassign[duplicated(tissue)==F][duplicated(cluster)==F][order(cluster)]
autoassign_v1

autoassign = rbind(autoassign_v1, autoassign_v0)[duplicated(cluster)==F]

# assignment_old = read_excel2(here("../1911_sympath_scRNA_r1lunge/data/zuordnung_SCT_snn_res.0.4_2020-01-21.xlsx"))
#
# stopifnot(all(lunge_merged$SCT_snn_res.0.4 %in% assignment$snn_res.0.4))
#
#
lunge_merged$tissue_v1 = autoassign[match_hk(lunge_merged$SCT_snn_res.0.4, autoassign$cluster), tissue]

table(lunge_merged$tissue_v1)

Idents(lunge_merged) %>% head

Idents(lunge_merged)  = lunge_merged$tissue_v1

p0 <- DimPlot(lunge_merged, reduction = "umap",pt.size = 0.8, split.by ="species", group.by = "SCT_snn_res.0.4", label = T, label.size = 5,repel = T, label.color = alpha("black", 0.5)) +NoLegend()
p0


jpeg(here("results/s507_4_three_species_clusternumbers_v2.jpg"), width = 24, height = 7.5, units = "in", res = 300, quality = 100)
print(p0)
dev.off()

p0b <- DimPlot(lunge_merged, reduction = "umap",pt.size = 0.9, group.by = "SCT_snn_res.0.4", label = T, label.size = 6, split.by = "species") +NoLegend()
p0b


p1 <- DimPlot(lunge_merged, reduction = "umap",pt.size = 0.7, group.by = "species")
p1


# p1 %>% ggplotly

p2 <- DimPlot(lunge_merged, reduction = "umap",pt.size = 1, group.by = "run10x")
p2


# p2 %>% ggplotly

p3 <- DimPlot(lunge_merged, reduction = "umap", ,pt.size = 1,group.by = "tissue")
p3


# p3 %>% ggplotly



p3b <- DimPlot(lunge_merged, reduction = "umap", group.by = "species_tissue")
p3b


# p3b  %>% ggplotly


p4 <- DimPlot(lunge_merged, reduction = "umap", group.by = "tissue", split.by = "species")
p4


# p4  %>%  ggplotly()


p5 =  DimPlot(lunge_merged, reduction = "umap",pt.size = 0.8, split.by ="species", group.by = "tissue", label = T, label.size = 3,repel = T, label.color = alpha("black", 1)) +NoLegend()
p5



jpeg(here("results/s507_4_three_species_originaltissueSUPPL.jpg"), width = 24, height = 7.5, units = "in", res = 300, quality = 100)
print(p5)
dev.off()
## png
##   2
lunge_merged$clustnum0.4_tissue_v1 = paste0(lunge_merged$SCT_snn_res.0.4, "_", lunge_merged$tissue_v1)
p5b =  DimPlot(lunge_merged, reduction = "umap",pt.size = 0.8, split.by ="species", group.by = "clustnum0.4_tissue_v1", label = T, label.size = 3,repel = T, label.color = alpha("black", 0.5)) +NoLegend()
p5b
## Warning: ggrepel: 2 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps


jpeg(here("results/s507_4_three_species_tissue_guessFrequency_v1.jpg"), width = 16, height = 5, units = "in", res = 300, quality = 100)
print(p5b)
dev.off()
## png
##   2
require(ggrepel)
head(p5b$data)

p5b_data_pre = data.table(p5b$data, keep.rownames = T)


p5b_data_pre[, median_x := median(UMAP_1), .(species, clustnum0.4_tissue_v1)]
p5b_data_pre[, median_y := median(UMAP_2), .(species, clustnum0.4_tissue_v1)]
p5b_data_label = unique(p5b_data_pre[,.(median_x,median_y, clustnum0.4_tissue_v1, species)])

p5c = ggplot(p5b_data_pre, aes(UMAP_1, UMAP_2, col = clustnum0.4_tissue_v1)) + geom_point(size = 1, alpha = 0.3) + facet_grid(.~species) + geom_text_repel(data = p5b_data_label, aes(x = median_x, y = median_y, label = clustnum0.4_tissue_v1), seed = 1902) + theme_cowplot() + guides(col = F) + geom_text_repel(data = p5b_data_label, aes(x = median_x, y = median_y, label = clustnum0.4_tissue_v1), seed = 1902, alpha = 0.4, col = "black") + Seurat:::FacetTheme()

p5c



jpeg(here("results/s507_4_three_species_tissue_guessFrequency_v1_Colored.jpg"), width = 24, height = 7.5, units = "in", res = 300, quality = 100)
print(p5c)
dev.off()

# speichern -----
lunge_merged@meta.data[1,]

lunge_merged$cluster_seurat_v1 = lunge_merged$SCT_snn_res.0.4
saveRDS(lunge_merged, file = here("results/s507_4_lunge_merged.RDS"))

# bb browser speichern----
sct_features = rownames(GetAssayData(lunge_merged[["SCT"]], slot = "counts"))
str(sct_features)
##  chr [1:33828] "1600002K03Rik" "1700020N01Rik" "1700109H08Rik" ...
lunge_mergedRNA = lunge_merged[sct_features,]


lunge_mergedRNA

lunge_mergedRNA <- SetAssayData(
  object = lunge_mergedRNA,
  slot = "counts",
  new.data = GetAssayData(lunge_merged[["SCT"]], slot = "counts"),
  assay = "RNA"
)

lunge_mergedRNA <- SetAssayData(
  object = lunge_mergedRNA,
  slot = "data",
  new.data = GetAssayData(lunge_merged[["SCT"]], slot = "data"),
  assay = "RNA"
)

lunge_mergedRNA <- SetAssayData(
  object = lunge_mergedRNA,
  slot = "scale.data",
  new.data = GetAssayData(lunge_merged[["SCT"]], slot = "scale.data"),
  assay = "RNA"
)

lunge_mergedRNA@assays$SCT = NULL
lunge_mergedRNA@assays

# Seurat::GetAssay(lunge_mergedRNA)
lunge_mergedRNA@active.assay

lunge_mergedRNA@active.assay = "RNA"
lunge_mergedRNA@active.assay

saveRDS(lunge_mergedRNA, file = here("results/s507_4_lunge_merged_NOsct_BBrowser.RDS"))



finalizeSkript()
