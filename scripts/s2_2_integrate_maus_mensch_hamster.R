#' ---
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'


# CREDITS: Parts of PIPELINE modified from https://github.com/immunogenomics/harmony/issues/41   kvshams commented on 11 Jun 2020 •

# # INITIATE SCRIPT ----
rm(list = ls())


r_on_server = T

if(r_on_server) {
  rootpath = "/net/ifs1/san_projekte/projekte/"
  .libPaths("/net/ifs1/san_projekte/projekte/genstat/07_programme/rpackages/forostar3")
  toolboxH::initializeSkript(computer = "forostar3" )
} else {
  rootpath =  "R:"
  toolboxH::initializeSkript(computer = "local" )
}
require(toolboxH) # https://github.com/holgerman/toolboxH
knitr::opts_chunk$set(cache = F, results = "markup", echo=T, fig.width = 10, fig.height = 7 )

packages2load =c("magrittr","knitr", "ggplot2", "scales", "ggthemes", "assertr","Seurat", 'paletteer',"pheatmap", "cowplot", "here", "plotly", "harmony")

new.packages <- packages2load[!(packages2load %in% installed.packages()[,"Package"])]
new.packages

if(length(new.packages)) BiocManager::install(new.packages)

for(i in packages2load) {

  suppressPackageStartupMessages(library(i, character.only = TRUE))

}

# # LOAD DATA ----

lung_list = readr::read_rds(here("results/s1_1_list_human_mouse_hamster_harmonized4integration.RDS"))


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


lunge_merged$species = factor(lunge_merged$species, levels = c("hamster", "mouse", "human charite", "human travaglini"))

lunge_merged$celltypeReported = lunge_merged$celltype
lunge_merged$celltype = NULL

lunge_merged$species_celltypeReported = paste0(lunge_merged$species,"_", lunge_merged$celltypeReported)

# # Some Plots----

p1 <- DimPlot(lunge_merged, reduction = "umap",pt.size = 0.7, group.by = "species")
p1

p2 <- DimPlot(lunge_merged, reduction = "umap",pt.size = 1, group.by = "run10x")
p2

# # SAVING----
saveRDS(lunge_merged, file = here("results/s2_1_lunge_merged.RDS"))
# lunge_merged = readr::read_rds(here("results/s2_1_lunge_merged.RDS"))


# # FINALIZE SCRIPT----
finalizeSkript()
