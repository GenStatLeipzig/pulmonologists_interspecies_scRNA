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

packages2load = c("magrittr","knitr", "ggplot2", "scales", "ggthemes", "assertr","Seurat", 'paletteer',"pheatmap", "cowplot", "here", "plotly")

new.packages <- packages2load[!(packages2load %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

for(i in packages2load) {

  suppressPackageStartupMessages(library(i, character.only = TRUE))

}

# # LOAD DATA ----
# Seurat objects from standard single species analysis, e.g. https://satijalab.org/seurat/articles/get_started.html,  QC filtered and set back applying Seurat::DietSeurat() function

human_travaglini  = readRDS(here("data/s511_1_diet_human_travaglini.RDS"))
human_travaglini

hamster = readRDS(here("data/s511_1_diet_hamster.RDS"))
hamster

mouse = readRDS(here("data/s511_1_diet_mouse.RDS"))
mouse

human_charite = readRDS(here("data/s511_1_diet_human_charite.RDS"))
human_charite



# # MAKE UNIQUE GENE NAMES for integration ----
gene_anno = fread(here("data/s512_1_renaming_table_orthologues_seurat.txt.gz"))

# downloaded from https://www.ensembl.org/biomart/martview/6dcbc1911b1c91a71d1d4055902e3c39  06/2021, Ensembl 104 GRCh38.p13
# reduced to genes existing in respective species specific seurat object and to unique assignments. In case of multiple assignments of one gene to another, i.e. multiple orthologues existed, the higher expressed gene was used


# # RENAME SEURAT OBJECTS using consensus name----
source(here("scripts/Seurat.update.gene.symbols.HGNC.R")) # from https://github.com/vertesy/Seurat.utils/

# ## rename hamster and check for proper renaming----
renaming_hamster = unique(gene_anno[species =="golden hamster", .(seurat_name, gene_used_in_consensus)])

renaming_hamster[allDuplicatedEntries(gene_used_in_consensus)] %>% assertr::verify(nrow(.)==0)
renaming_hamster[allDuplicatedEntries(seurat_name)] %>% assertr::verify(nrow(.)==0)
showNA(renaming_hamster)%>% assertr::verify(sum(.$NAs)==0)
stopifnot(all(renaming_hamster$seurat_name %in% rownames(hamster)))

hamster3=hamster[renaming_hamster$seurat_name,]
hamster3
grep("^NA$", rownames(hamster3))%>% assertr::verify(length(.)==0)

stopifnot(identical(dimnames(hamster3@assays$RNA)[[1]], renaming_hamster$seurat_name))

hamster3 = RenameGenesSeurat(obj = hamster3, newnames = renaming_hamster$gene_used_in_consensus)
hamster3
grep("^NA$", rownames(hamster3))%>% assertr::verify(length(.)==0)
rownames(hamster3)[duplicated(rownames(hamster3))]%>% assertr::verify(length(.)==0)


# ## rename mouse and check for proper renaming----
renaming_mouse = unique(gene_anno[species =="mouse", .(seurat_name, gene_used_in_consensus)])

renaming_mouse[allDuplicatedEntries(seurat_name)] %>% assertr::verify(nrow(.)==0)
renaming_mouse[allDuplicatedEntries(gene_used_in_consensus)] %>% assertr::verify(nrow(.)==0)
showNA(renaming_mouse)%>% assertr::verify(sum(.$NAs)==0)
stopifnot(all(renaming_mouse$seurat_name %in% rownames(mouse)))

mouse3=mouse[renaming_mouse$seurat_name,]
mouse3
grep("^NA$", rownames(mouse3))%>% assertr::verify(length(.)==0)
stopifnot(identical(dimnames(mouse3@assays$RNA)[[1]], renaming_mouse$seurat_name))

mouse3 = RenameGenesSeurat(obj = mouse3, newnames = renaming_mouse$gene_used_in_consensus)
mouse3
grep("^NA$", rownames(mouse3))%>% assertr::verify(length(.)==0)
rownames(mouse3)[duplicated(rownames(mouse3))] %>% assertr::verify(length(.)==0)


# # ADD SPECIES, BATCH, AND CELLTYPE consistently ----
# in each Seurat object Seurats SCT workflow must be done in each batch separately
human_travaglini @meta.data$run10x = paste0("human_",human_travaglini @meta.data$channel)
table(human_travaglini @meta.data$run10x)

human_travaglini @meta.data$species = "human travaglini"

human_charite@meta.data$run10x = paste0("human_",human_charite@meta.data$orig.ident)
table(human_charite@meta.data$run10x)

human_charite@meta.data$species = "human charite"


mouse3@meta.data$run10x = paste0("mouse_",mouse3@meta.data$orig.ident)
table(mouse3@meta.data$run10x)

mouse3@meta.data$species = "mouse"

hamster3@meta.data$run10x = paste0("hamster_",hamster3@meta.data$orig.ident)
table(hamster3@meta.data$run10x)

hamster3@meta.data$species = "hamster"

# ## Add the originally reported celltype  ----

human_travaglini$celltype = human_travaglini$free_annotation
human_charite$celltype =  human_charite$predicted.id
table(mouse3$celltype)
table(hamster3$celltype)



# # CREATE A SINGLE scRNA-SEQ list for later use----

human_travaglini_list <- SplitObject(object = human_travaglini , split.by = "run10x")
human_charite_list <- SplitObject(object = human_charite, split.by = "run10x")
hamster_list <- SplitObject(object = hamster3, split.by = "run10x")
lung_list = c(hamster_list, mouse3, human_travaglini_list, human_charite_list)
names(lung_list)
names(lung_list)[4] = unique(mouse3$run10x)
names(lung_list)

rm(hamster_list)
rm(human_travaglini_list)
rm(human_charite_list)


# # REMOVE - IF EXISTENT - NON RNA assays----
for (i in 1:length(lung_list)) {
  message("Removing other assays than RNA from entry ", i, " i.e. ", names(lung_list)[i], "\n-------------------------------------------------------")
  mydata=lung_list[[i]]
  try(mydata[['integrated']] <- NULL,silent = T)
  try(mydata[['SCT']] <- NULL, silent = T)


  lung_list[[i]] <- mydata
}


# # STATS----
dim(human_travaglini )
dim(human_charite)

dim(mouse)
dim(mouse3)
dim(mouse)[1] - dim(mouse3)[1]
dim(mouse3)[1]/dim(mouse)[1]

dim(hamster)
dim(hamster3)
dim(hamster)[1] - dim(hamster3)[1]
dim(hamster3)[1]/dim(hamster)[1]


# # SCT TRANSFORM ----
for (i in 1:length(lung_list)) {
  message("SCT on ", i, " i.e. ", names(lung_list)[i], "\n-------------------------------------------------------")

  lung_list[[i]] <- SCTransform(lung_list[[i]], verbose = T)
}


print(table(names(warnings() ))) # "iteration limit reached" see comment https://github.com/ChristophH/sctransform/issues/25 :"These warnings are showing that there are some genes for which it is hard to reliably estimate theta (presumably because of very few non-zero observations). Usually we don't worry about these warnings too much, since we regularize the parameters in a later step, thus averaging out uncertainty of individual gene  parameters."

# # CREATE SINGLE scRNA-Seq OBJECT ----


lunge_merged <- merge(lung_list[[1]],
                      y = lung_list[2:length(lung_list)],
                      project = "lung",
                      merge.data = TRUE)

if(DefaultAssay(lunge_merged) != "SCT") DefaultAssay(lunge_merged) = "SCT"


# # PCA-DIMENSION REDUCTION----
list.features <- SelectIntegrationFeatures(object.list = lung_list, nfeatures = 3000)
VariableFeatures(lunge_merged) <- list.features

lunge_merged <- RunPCA(object = lunge_merged, assay = "SCT", npcs = 50 ,features = list.features)


# # BATCH REDUCTION via HARMONY-----

lunge_merged <- RunHarmony(object = lunge_merged,
                           assay.use = "SCT",
                           reduction = "pca",
                           dims.use = 1:50,
                           group.by.vars = "run10x",
                           plot_convergence = TRUE)

# # UMAP AND CLUSTERING -----

lunge_merged <- RunUMAP(object = lunge_merged, assay = "SCT", reduction = "harmony", dims = 1:50)


lunge_merged <- FindNeighbors(object = lunge_merged, assay = "SCT", reduction = "harmony", dims = 1:50)
#
lunge_merged <- FindClusters(object = lunge_merged, resolution = 0.4)
table(lunge_merged$SCT_snn_res.0.4)

dput(unique(lunge_merged$species))


lunge_merged$species = factor(lunge_merged$species, levels = c( "mouse", "hamster","human charite", "human travaglini"))

lunge_merged$celltypeReported = lunge_merged$celltype

relevantnames = c("orig.ident", "nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT",
                  "celltypeReported", "run10x", "species", "SCT_snn_res.0.4")

lunge_merged@meta.data = lunge_merged@meta.data[,relevantnames]

# # Rename Clusters----
# using original assignments with some marker gene based corrections


celltypeassignment = c("0" = "T Cells", '1' = "Endothelial", '2' = "Alveolar Macrophages", '3' = "NK Cells", '4' = "Macrophages/ Monocytes", '5' = "Endothelial", '6' = "Fibroblasts", '7' = "AT2", '8' = "Perivascular", '9' = "B Cells", '10' = "Macrophages/ Monocytes", '11' = "Mast Cells", '12' = "Macrophages/ Monocytes", '13' = "Endothelial", '14' = "AT1", '15' = "Alveolar Macrophages", '16' = "Ciliated Cells", '17' = "Endothelial", '18' = "Club Cells", "19" = "ly Endothelial", "20" = "Dendritic Cells", "21" = "Endothelial", "22" = "T Cells", "23" = "Proliferating NK/T", "24" = "Proliferating Alveolar Macrophages", "25" = "Endothelial", "26" = "T Cells", "27" = "B Cells")


clusternames_new = data.table(old_id  = names(celltypeassignment),
                              new_id = celltypeassignment)

lunge_merged$cluster_seurat_v2 = clusternames_new[match_hk(lunge_merged$SCT_snn_res.0.4 %>% as.character, clusternames_new$old_id%>% as.character), new_id]


clusterorder = c('Alveolar Macrophages', 'Proliferating Alveolar Macrophages', 'Macrophages/ Monocytes', 'Mast Cells', 'Dendritic Cells', 'T Cells', 'NK Cells', 'Proliferating NK/T', 'B Cells', 'AT1', 'AT2', 'Ciliated Cells', 'Club Cells', 'Endothelial', 'ly Endothelial', 'Fibroblasts', 'Perivascular' )

lunge_merged$cluster_seurat_v2 = factor(lunge_merged$cluster_seurat_v2, levels = clusterorder)


# # SAVING----
saveRDS(lunge_merged, file = here("results/s1_3_lunge_merged.RDS"))



# # FINALIZE SCRIPT----
finalizeSkript()


