rm(list = ls())
require(toolboxH)
require(here)
require(Seurat)
require(ggplot2)
require(ggthemes)
require(scales)
require(dplyr)
require(tidyr)
require(patchwork)

require(harmony)

# Load gene annotation accross species
orthologues = fread(here("results/s602_9v2_orthologues.txt.gz")) # "Ensembl release 104 - May 2021 Human genes (GRCh38.p13), https://www.ensembl.org/biomart/martview/bf27d7b2470fe6d55826530cdb586229 processed in previous script
orthologues

mouse = readRDS(here("results/s602_9v2_mouse_preprocessed_withDoubletts.RDS"))
mouse = DietSeurat(mouse)

monkey = readRDS(here("results/s602_9v2_monkey_preprocessed_withDoubletts.RDS"))
monkey = DietSeurat(monkey)

pig = readRDS(here("results/s602_9v2_pig_preprocessed_withDoubletts.RDS"))
pig = DietSeurat(pig)

rat = readRDS(here("results/s602_9v2_rat_preprocessed_withDoubletts.RDS"))
rat = DietSeurat(rat)

hamster = readRDS(here("results/s602_9v2_hamster_preprocessed_withDoubletts.RDS"))
hamster = DietSeurat(hamster)

humanTravaglini = readRDS(here("results/s602_9v2_humanTravaglini10x_preprocessed_withDoubletts.RDS"))
humanTravaglini =   DietSeurat(humanTravaglini)

humanCharite = readRDS(here("results/s602_9v2_humanCharite_preprocessed_withDoubletts.RDS"))
humanCharite =   DietSeurat(humanCharite)

humanSS2 = readRDS(here("results/s602_9v2_humanTravagliniSS2_preprocessed_withDoubletts.RDS"))
humanSS2 =   DietSeurat(humanSS2)


# >> harmonizing gene names ----
# Ensure that  "" is replaced by  NA in Gene names----
showNA(orthologues, showAllNoNA = F)
for(i in c("Pig gene name", "Rat gene name", "Mouse gene name", "Golden Hamster gene name",
           "Vervet-AGM gene name", "Human gene name", "Rat gene name uc")) {
  orthologues[get(i)=="",(i):=NA]
}
showNA(orthologues, showAllNoNA = F)


# use either name or if not present stable ID ----

orthologues[,name_human := ifelse(is.na(`Human gene name`), `Human gene stable ID`, `Human gene name`)]

orthologues[, name_hamster := ifelse(`Golden Hamster gene name` %nin% rownames(hamster@assays$RNA),`Golden Hamster gene stable ID`, `Golden Hamster gene name` )]

orthologues[, name_monkey := ifelse(`Vervet-AGM gene name` %nin% rownames(monkey@assays$RNA),`Vervet-AGM gene stable ID`,`Vervet-AGM gene name`)]

orthologues[, name_pig := ifelse(`Pig gene name` %nin% rownames(pig@assays$RNA),`Pig gene stable ID`, `Pig gene name` )]


orthologues[, name_rat := ifelse(`Rat gene name uc` %nin% rownames(rat@assays$RNA),`Rat gene stable ID`, `Rat gene name uc` )]
table(orthologues$name_rat %in% unique(rownames(rat@assays$RNA)))

orthologues[, name_mouse := ifelse(`Mouse gene name` %nin% rownames(mouse@assays$RNA),`Mouse gene stable ID`, `Mouse gene name` )]


# count frequency of available orthologues in database ----
orthologues[,n_genenames := as.numeric(name_human !="") +
              as.numeric(name_mouse !="") +
              as.numeric(name_hamster !="") +
              as.numeric(name_monkey !="") +
              as.numeric(name_pig !="")+
              as.numeric(name_rat !="")


]
count_orthologue =  orthologues[,uniqueN(`Human gene stable ID`, na.rm = T),n_genenames]
count_orthologue

message("In Ensembl release 104 - May 2021 Human genes (GRCh38.p13), a total of ", count_orthologue[n_genenames==6, V1], " genes have orthologues in all six species.")

# Ensure that replacement with Gene ID did not introduce "" instead of NA----
showNA(orthologues, showAllNoNA = F)
for(i in c("name_hamster",  "name_human", "name_mouse","name_monkey","name_pig", "name_rat")) {
  orthologues[get(i)=="",(i):=NA]
}
showNA(orthologues, showAllNoNA = F)

# prioritize for duplicated genes higher expressed gene ----

av_express_humanCharite = AverageExpression(humanCharite,assays = "RNA", slot = "count")
av_express_humanCharite_dt = data.table(name_human = rownames(av_express_humanCharite$RNA),
                                         max_in_any_mensch = apply(av_express_humanCharite$RNA, 1, max))

av_express_humanCharite_dt

av_express_humanTravaglini = AverageExpression(humanTravaglini,assays = "RNA", slot = "count")
av_express_humanTravaglini_dt = data.table(name_human = rownames(av_express_humanTravaglini$RNA),
                                            max_in_any_mensch = apply(av_express_humanTravaglini$RNA, 1, max))

av_express_humanTravaglini_dt


av_express_humanSS2 = AverageExpression(humanSS2,assays = "RNA", slot = "count")
av_express_humanSS2_dt = data.table(name_human = rownames(av_express_humanSS2$RNA),
                                           max_in_any_mensch = apply(av_express_humanSS2$RNA, 1, max))

av_express_humanSS2_dt



orthologues[, gx_humanCharite := av_express_humanCharite_dt[match_hk(orthologues$name_human, av_express_humanCharite_dt$name_human),max_in_any_mensch]]

orthologues[, gx_humanTravaglini := av_express_humanTravaglini_dt[match_hk(orthologues$name_human, av_express_humanTravaglini_dt$name_human),max_in_any_mensch]]

orthologues[, gx_humanSS2 := av_express_humanSS2_dt[match_hk(orthologues$name_human, av_express_humanSS2_dt$name_human),max_in_any_mensch]]

orthologues[,max_in_any_mensch := ifelse(is.na(gx_humanCharite) & is.na(gx_humanTravaglini), gx_humanSS2,
                                         ifelse(is.na(gx_humanCharite), gx_humanTravaglini, gx_humanCharite))]


# reating table for renaming gene names ----
require(Seurat.utils) #  # from https://github.com/vertesy/Seurat.utils/

renameSeuratOrthologues = function(seuratobject, speciesname, orthotable,speciescolumn, confidencecolumn , all_genes_uppercase =F) {
  # seuratobject = hamster; speciesname = "hamster"; orthotable = copy(orthologues); confidencecolumn = "Golden Hamster orthology confidence [0 low, 1 high]"; speciescolumn = "name_hamster"
  #
  orthologues_species= unique(orthotable[is.na(get(speciescolumn))==F, .(name_human, name_species = get(speciescolumn), max_in_any_mensch,confidence=get(confidencecolumn))])[order(name_species, -max_in_any_mensch, na.last = T)]
  orthologues_species[allDuplicatedEntries(name_species)]
  orthologues_species= orthologues_species[duplicated(name_species)==F]
  orthologues_species

  av_express_species= AverageExpression(seuratobject,assays = "RNA", slot = "count")
  hh(av_express_species$RNA)
  dim(av_express_species$RNA)
  renaming_species= data.table(species = speciesname, name_original = rownames(av_express_species$RNA),
                               max_in_any_celltype = apply(av_express_species$RNA, 1, max))

  if(all_genes_uppercase==T) {
    renaming_species[, name_human := orthologues_species[match_hk(renaming_species$name_original %>% toupper(), orthologues_species$name_species %>% toupper()), name_human]]
    renaming_species[, confidence := orthologues_species[match_hk(renaming_species$name_original %>% toupper(), orthologues_species$name_species %>% toupper()), confidence]]
  } else {
    renaming_species[, name_human := orthologues_species[match_hk(renaming_species$name_original, orthologues_species$name_species), name_human]]
    renaming_species[, confidence := orthologues_species[match_hk(renaming_species$name_original, orthologues_species$name_species), confidence]]
  }

  renaming_species= renaming_species[, gene_used_in_consensus := ifelse(is.na(name_human), name_original, name_human)][order(-name_human, -max_in_any_celltype, na.last = T)]
  renaming_species[allDuplicatedEntries(gene_used_in_consensus)]
  renaming_species[, genes2remove := duplicated(name_human)]
  renaming_species[, .N,genes2remove]
  message("Mean confidence ", speciesname, ":",renaming_species[genes2remove==F & is.na(name_human)==F, mean(confidence)])
  print(renaming_species[genes2remove==F & is.na(name_human)==F, table(confidence)])

  seuratobject = seuratobject[renaming_species[genes2remove==F, name_original],]
  seuratobject
  grep("^NA$", rownames(seuratobject))%>% assertr::verify(length(.)==0)
  stopifnot(identical(dimnames(seuratobject@assays$RNA)[[1]], renaming_species[genes2remove==F, name_original ]))


  seuratobject = RenameGenesSeurat(obj = seuratobject, newnames = renaming_species[genes2remove==F, gene_used_in_consensus])
  grep("^NA$", rownames(seuratobject))%>% assertr::verify(length(.)==0)
  res = c()
  res$seuratobject = seuratobject
  res$renamingtable = renaming_species
  res
}


hamster_renamelist = renameSeuratOrthologues(seuratobject = hamster,
                                             speciesname = "hamster",
                                             orthotable = orthologues,
                                             speciescolumn = "name_hamster",
                                             confidencecolumn = "Golden Hamster orthology confidence [0 low, 1 high]",
                                             all_genes_uppercase =F)


mouse_renamelist = renameSeuratOrthologues(seuratobject = mouse,
                                           speciesname = "mouse",
                                           orthotable = orthologues,
                                           speciescolumn = "name_mouse",
                                           confidencecolumn = "Mouse orthology confidence [0 low, 1 high]",
                                           all_genes_uppercase =F)

monkey_renamelist = renameSeuratOrthologues(seuratobject = monkey,
                                            speciesname = "monkey",
                                            orthotable = orthologues,
                                            speciescolumn = "name_monkey",
                                            confidencecolumn = "Vervet-AGM orthology confidence [0 low, 1 high]",
                                            all_genes_uppercase =F)


pig_renamelist = renameSeuratOrthologues(seuratobject = pig,
                                         speciesname = "pig",
                                         orthotable = orthologues,
                                         speciescolumn = "name_pig",
                                         confidencecolumn = "Pig orthology confidence [0 low, 1 high]",
                                         all_genes_uppercase =F)




rat_renamelist = renameSeuratOrthologues(seuratobject = rat,
                                         speciesname = "rat",
                                         orthotable = orthologues,
                                         speciescolumn = "name_rat",
                                         confidencecolumn = "Rat orthology confidence [0 low, 1 high]",
                                         all_genes_uppercase =T)





# ##count realized gene overlap ----

makeInputUpset = function(plotlist) {
  library(data.table)
  if(length(names(plotlist))==0) {
    message("no names for list entries found -  providing standardnames  `list1...")
    names_list = paste0("list", seq(along = plotlist))
    names(plotlist) = names_list
  } else names_list = names(plotlist)

  plotlist2 = lapply(names_list, function(myname) {
    data.table(variable = myname,value = plotlist[[myname]])
  }
  )

  plotlist3 = rbindlist(plotlist2)

  plotlist4 = dcast.data.table(plotlist3, value ~ variable, fun.aggregate = function(x) as.numeric(length(x)>0 ))
  plotlist4

}


require(UpSetR)


input1 = makeInputUpset(list(hamster = hamster_renamelist$renamingtable$gene_used_in_consensus,
                             humanCharite = rownames(humanCharite),
                             humanTravaglini = rownames(humanTravaglini),
                             humanSS2 = rownames(humanSS2),
                             monkey = monkey_renamelist$renamingtable$gene_used_in_consensus,
                             mouse = mouse_renamelist$renamingtable$gene_used_in_consensus,
                             pig = pig_renamelist$renamingtable$gene_used_in_consensus,
                             rat = rat_renamelist$renamingtable$gene_used_in_consensus))


upset(input1,nsets = ncol(input1)-1,order.by =c("freq" ,"degree"), decreasing = c(TRUE, TRUE), set_size.scale_max =60000, text.scale = 1.3, set_size.show = TRUE,mainbar.y.label = "Overlap Genes", nintersects = 8)

jpeg(here("results/s603_2_overlap_genes_all_datasets.jpeg"), width = 6.5, height = 4.5, units = "in", res = 150, quality = 100)

upset(input1,nsets = ncol(input1)-1,order.by =c("freq" ,"degree"), decreasing = c(TRUE, TRUE), set_size.scale_max =60000, text.scale = 1.3, set_size.show = TRUE,mainbar.y.label = "Overlap Genes", nintersects = 8)

dev.off()

pdf(here("results/s603_2_overlap_genes_all_datasets.pdf"), width = 6.5, height = 4.5)

upset(input1,nsets = ncol(input1)-1,order.by =c("freq" ,"degree"), decreasing = c(TRUE, TRUE), set_size.scale_max =60000, text.scale = 1.3, set_size.show = TRUE,mainbar.y.label = "Overlap Genes", nintersects = 8)

dev.off()

# high confidence overlap
input3 = makeInputUpset(list(hamster = hamster_renamelist$renamingtable[confidence ==1, gene_used_in_consensus],
                             humanCharite = rownames(humanCharite),
                             humanTravaglini = rownames(humanTravaglini),
                             humanSS2 = rownames(humanSS2),
                             monkey = monkey_renamelist$renamingtable[confidence ==1, gene_used_in_consensus],
                             mouse = mouse_renamelist$renamingtable[confidence ==1, gene_used_in_consensus],
                             pig = pig_renamelist$renamingtable[confidence ==1, gene_used_in_consensus],
                             rat = rat_renamelist$renamingtable[confidence ==1, gene_used_in_consensus]))


upset(input3,nsets = ncol(input3)-1,order.by =c("freq" ,"degree"), decreasing = c(TRUE, TRUE), set_size.scale_max =60000, text.scale = 1.3, set_size.show = TRUE,mainbar.y.label = "Overlap Genes\n(high orthologue confidence)", nintersects = 8)

renamelist_all = rbind(hamster_renamelist$renamingtable,
                       mouse_renamelist$renamingtable,
                       monkey_renamelist$renamingtable,
                       pig_renamelist$renamingtable,
                       rat_renamelist$renamingtable)

renamelist_all

## add stressgenes ----
stressgene_pre = fread(here("data/stress_genes.txt"))
stressgene_pre

stressgene_pre2 = mouse_renamelist$renamingtable[name_original %in%stressgene_pre$genename_harmonized]
stressgene_pre2

stressgene = renamelist_all[genes2remove ==F &name_human %in% stressgene_pre2$name_human ]
stressgene2 = dcast.data.table(stressgene, name_human ~ species, value.var = "name_original")
stressgene2
stressgene_complete = stressgene2[is.na(pig)==F &
                                    is.na(hamster)==F &
                                    is.na(monkey)==F &
                                    is.na(mouse )==F &
                                    is.na(rat)==F]

all(stressgene_complete$name_human %in% rownames(humanCharite))
all(stressgene_complete$name_human %in% rownames(humanTravaglini))

stressgene_complete$name_human%>% sort %>% paste(., collapse = ", ")

qlit341 = venn4(stressgene_complete$name_human, rownames(humanCharite), rownames(humanTravaglini), rownames(humanSS2)) # if stressgenes are complete in species, they are also found in our three human atasets

stressgene2[,complete := name_human%in% stressgene_complete$name_human]

hamster2 = hamster_renamelist$seuratobject
mouse2 = mouse_renamelist$seuratobject
monkey2 = monkey_renamelist$seuratobject
pig2 = pig_renamelist$seuratobject
rat2 = rat_renamelist$seuratobject

intersectgenes = Reduce(intersect, list(hamster = hamster_renamelist$renamingtable$gene_used_in_consensus,
                                        humanCharite = rownames(humanCharite),
                                        humanTravaglini = rownames(humanTravaglini),
                                        humanSS2 = rownames(humanSS2),
                                        monkey = monkey_renamelist$renamingtable$gene_used_in_consensus,
                                        mouse = mouse_renamelist$renamingtable$gene_used_in_consensus,
                                        pig = pig_renamelist$renamingtable$gene_used_in_consensus,
                                        rat = rat_renamelist$renamingtable$gene_used_in_consensus)) # https://www.r-bloggers.com/2012/06/intersect-for-multiple-vectors-in-r/

str(intersectgenes)

fwrite(orthologues, here("results/s603_2_orthologues_used.txt.gz"), sep = "\t")
fwrite(renamelist_all, here("results/s603_2_renamelist_all.txt.gz"), sep = "\t")


fwrite( stressgene2, here("results/s603_2_stressgenes_in__species.txt"))

# save(hamster2, mouse2, monkey2, pig2, rat2, humanCharite, humanTravaglini, intersectgenes, file = here("results/s603_2_renamed_input_seurat.RData"))
#
#
# finalizeSkript()
# finalizeSkript()



# ## renaming  the originally reported celltype for easier comparison  ----

humanTravaglini$celltype = humanTravaglini$free_annotation
humanCharite$celltype =  humanCharite$predicted.id
humanSS2$celltype = humanSS2$free_annotation
table(mouse2$celltype)
table(hamster2$celltype)

monkey2$celltype = "not available"
pig2$celltype = "not available"
rat2$celltype = "not available"


# # CREATE A SINGLE scRNA-SEQ list for later use----

humanTravaglini$pct.mito = mean(humanSS2$pct.mito)
humanTravaglini_list <- SplitObject(object = humanTravaglini , split.by = "run10x")

humanCharite_list <- SplitObject(object = humanCharite, split.by = "run10x")

hamster_list <- SplitObject(object = hamster2, split.by = "run10x")

monkey_list = SplitObject(object = monkey2, split.by = "run10x")

pig_list = SplitObject(object = pig2, split.by = "run10x")

rat_list = SplitObject(object = rat2, split.by = "run10x")

lung_list = c(hamster_list, mouse2, humanTravaglini_list, humanCharite_list, monkey_list, pig_list, rat_list, humanSS2)
names(lung_list)
names(lung_list)[4] = unique(mouse2$run10x)
names(lung_list)[23] = unique(humanSS2$run10x)
names(lung_list)

# # REMOVE - IF EXISTENT - NON RNA assays----
for (i in 1:length(lung_list)) {
  message("Removing other assays than RNA from entry ", i, " i.e. ", names(lung_list)[i], "\n-------------------------------------------------------")
  mydata=lung_list[[i]]
  try(mydata[['integrated']] <- NULL,silent = T)
  try(mydata[['SCT']] <- NULL, silent = T)


  lung_list[[i]] <- mydata
}


# here's some basic statistics before final filtering
sobj <- merge(lung_list[[1]],
              y = lung_list[2:length(lung_list)],
              project = "lung",
              merge.data = TRUE)


sobj[['log10nCount_RNA']] <- log10(sobj[['nCount_RNA']])
stats1  = sobj@meta.data %>%
  dplyr::group_by(run10x) %>%
  dplyr::summarise(ncells=dplyr::n(),
                   ndoublets=sum(DF.classifications_ConsidHomoDoubl=="Singlet"),
                   nGene=median(nFeature_RNA),
                   nUMI=median(nCount_RNA),
                   pct.mito=mean(pct.mito),
                   pct.ribo=mean(pct.ribo))
stats1

ggplot(sobj@meta.data %>%
         dplyr::select(run10x, log10nCount_RNA, nFeature_RNA,  pct.mito,
                       pct.ribo) %>%
         gather(metric,value,-run10x),
       aes(x=run10x,y=value,fill=run10x)) +
  geom_boxplot(outlier.size=.5) +
  facet_wrap(~metric,ncol=4,scales='free_y') +
  theme(axis.text.x=element_blank())


# # >> REmove Doublettes and ensure basic QC ----

param_QC_filter_mito_max 	=30 # ABX 30, human 10
param_QC_filter_RNAfeatures_min 	=300
param_QC_filter_RNAfeatures_max 	= 6000
param_QC_filter_RNAcounts_min 	=1000
param_QC_filter_RNAcounts_max	=35000

param_QC_filter_RNAcountsSS2_max	=4000000





for (i in 1:length(lung_list)) {
  message("Final QC and Doublett filtering of experiment", i, " i.e. ", names(lung_list)[i], "\n-------------------------------------------------------")
  mydata=lung_list[[i]]
  dim_before = dim(mydata)
  myname = names(lung_list[i])

  if(grepl("SS2", myname)==F) {
        mydata = subset(mydata, (pct.mito < param_QC_filter_mito_max) &
           (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
           (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
           (nCount_RNA < param_QC_filter_RNAcounts_max) &
           (nCount_RNA > param_QC_filter_RNAcounts_min) &
           DF.classifications_ConsidHomoDoubl == "Singlet")
  } else {
    mydata = subset(mydata, (pct.mito < param_QC_filter_mito_max) &
                      (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                      (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                      (nCount_RNA < param_QC_filter_RNAcountsSS2_max) &
                      (nCount_RNA > param_QC_filter_RNAcounts_min) &
                      DF.classifications_ConsidHomoDoubl == "Singlet")



  }
  dim_after = dim(mydata)
  message(Sys.time(),"...Filtering ", myname, " from\n", dim_before[1], " genes x ", dim_before[2], " cells to \n", dim_after[1], " genes x ", dim_after[2], " cells")
  message("(i.e. removing ", ((dim_before[2]-dim_after[2])/dim_after[2]) %>% proz()," of given cells)")

  lung_list[[i]] <- mydata
}

# here's some basic statistics before final filtering
sobj <- merge(lung_list[[1]],
              y = lung_list[2:length(lung_list)],
              project = "lung",
              merge.data = TRUE)


sobj[['log10nCount_RNA']] <- log10(sobj[['nCount_RNA']])
stats2  = sobj@meta.data %>%
  dplyr::group_by(run10x) %>%
  dplyr::summarise(ncells=dplyr::n(),
                   ndoublets=sum(DF.classifications_ConsidHomoDoubl=="Singlet"),
                   nGene=median(nFeature_RNA),
                   nUMI=median(nCount_RNA),
                   pct.mito=mean(pct.mito),
                   pct.ribo=mean(pct.ribo))
stats2

qcplot2 =ggplot(sobj@meta.data %>%
         dplyr::select(run10x, log10nCount_RNA, nFeature_RNA,  pct.mito,
                       pct.ribo) %>%
         gather(metric,value,-run10x),
       aes(x=run10x,y=value,fill=run10x)) +
  geom_boxplot(outlier.size=.5) +
  facet_wrap(~metric,ncol=4,scales='free_y') +
  theme(axis.text.x=element_blank())

qcplot2

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
require(harmony)
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

lunge_merged <- FindClusters(object = lunge_merged, resolution = 0.3)
table(lunge_merged$SCT_snn_res.0.3)

dput(unique(lunge_merged$species))


table(lunge_merged$run10x, lunge_merged$species,useNA = "always")

lunge_merged$species = factor(lunge_merged$species, levels = c("monkey", "pig", "rat", "mouse", "hamster","humanCharite", "humanTravaglini", "humanSS2"))

lunge_merged$celltypeReported = lunge_merged$celltype



# plotten ----
DimPlot(lunge_merged, group.by = "SCT_snn_res.0.4")
DimPlot(lunge_merged, group.by = "SCT_snn_res.0.3")
DimPlot(lunge_merged, group.by = "SCT_snn_res.0.4", split.by = "species", ncol = 4)
DimPlot(lunge_merged, group.by = "run10x")



jpeg(here("results/s603_2_lunge_merged_SCT_snn_res.0.4_all.jpeg"), width = 7,7, res= 150, quality = 100, units = "in")
DimPlot(lunge_merged, group.by = "SCT_snn_res.0.4")
dev.off()


jpeg(here("results/s603_2_lunge_merged_run10x_all.jpeg"), width = 7,7, res= 150, quality = 100, units = "in")
DimPlot(lunge_merged, group.by = "run10x")
dev.off()


jpeg(here("results/s603_2_lunge_merged_SCT_snn_res.0.3_all.jpeg"), 7,7, res= 150, quality = 100, units = "in")
DimPlot(lunge_merged, group.by = "SCT_snn_res.0.3")
dev.off()

jpeg(here("results/s603_2_lunge_merged_SCT_snn_res.0.4_by_species.jpeg"), 14,7, res= 150, quality = 100, units = "in")
DimPlot(lunge_merged, group.by = "SCT_snn_res.0.4", split.by = "species", ncol = 4)
dev.off()

jpeg(here("results/s603_2_lunge_merged_SCT_snn_res.0.3_by_species.jpeg"), 14,7, res= 150, quality = 100, units = "in")
DimPlot(lunge_merged, group.by = "SCT_snn_res.0.3", split.by = "species", ncol = 4)
dev.off()


FeaturePlot(lunge_merged, c('EPCAM','COL1A2','CLDN5','PTPRC',
                    'KRT5','MUC5B'),
            ncol=3, sort=TRUE, label=TRUE)

run_markersearch = T
if(run_markersearch==T) { # runs long
markers <- FindAllMarkers(lunge_merged, only.pos = TRUE, min.pct = 0.25,
                          thresh.use = 0.25,
                          print.bar=FALSE, verbose=FALSE)

WriteXLS_hk('markers',here('results/s603_2_lung_integrated_markers.xlsx'))


top5 <- markers %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC )

top10 <- markers %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC )

DotPlot(lunge_merged, features = unique(top5$gene)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
}

run_markers_res0_3_search =T
if(run_markers_res0_3_search==T) { # runs long

  lunge_merged <- SetIdent(lunge_merged,  value = 'SCT_snn_res.0.3')

  markers_res0_3 <- FindAllMarkers(lunge_merged, only.pos = TRUE, min.pct = 0.25,
                          thresh.use = 0.25,
                          print.bar=FALSE, verbose=FALSE)

WriteXLS_hk('markers_res0_3',here('results/s603_2_lung_integrated_markers.xlsx'))


top5 <- markers_res0_3 %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC )

top10 <- markers_res0_3 %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC )

DotPlot(lunge_merged, features = unique(top5$gene)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
}

# # Compare  Clusters with previous object----


lung_list_1st_subm = readRDS(here("data/s1_3_lunge_merged.RDS"))
Seurat::Assays(lung_list_1st_subm)
DimPlot(lung_list_1st_subm)


lunge_merged$old_celltype = lung_list_1st_subm@meta.data[match_hk(colnames(lunge_merged),colnames(lung_list_1st_subm)),"cluster_seurat_v2"]

DimPlot(lunge_merged, group.by = "old_celltype")

helper = lunge_merged@meta.data %>% as.data.table(keep.rownames = T)
helper2 = helper[,.N, .(SCT_snn_res.0.4,old_celltype)]
helper3 = helper2[is.na(old_celltype)==F][order(-N)][duplicated(SCT_snn_res.0.4)==F]
venn2(helper3$SCT_snn_res.0.4, lunge_merged$SCT_snn_res.0.4)

relevantnames = c("orig.ident", "nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT",'pct.ribo' , 'pct.mito',"celltypeReported", "run10x", "species", "SCT_snn_res.0.3", "SCT_snn_res.0.4", 'old_celltype')

lunge_merged@meta.data = lunge_merged@meta.data[,relevantnames]

venn2(relevantnames,names(lunge_merged@meta.data))

toolboxH::plotSankey(lunge_merged@meta.data[,c('SCT_snn_res.0.3', 'SCT_snn_res.0.4')], spalte4color = "SCT_snn_res.0.4")


lunge_merged$dataset = lunge_merged$species
lunge_merged$species = ifelse(grepl("human", lunge_merged$dataset), "human", lunge_merged$dataset %>% as.character())
lunge_merged$species %>% table()


table(lunge_merged$SCT_snn_res.0.4)
lunge_merged@meta.data$old_celltype_completed = helper3[match_hk(lunge_merged$SCT_snn_res.0.4, helper3$SCT_snn_res.0.4),old_celltype %>% as.character()]

lunge_merged@meta.data$old_celltype_completed = ifelse(is.na(lunge_merged$old_celltype_completed), lunge_merged$SCT_snn_res.0.4 %>% as.character(), lunge_merged$old_celltype_completed)
mytable(lunge_merged@meta.data$old_celltype_completed)

DimPlot(lunge_merged, group.by = "old_celltype_completed", label = T, repel = T) + guides(color = "none")
DimPlot(lunge_merged, group.by = "old_celltype_completed") %>% plotly::ggplotly()

DimPlot(lunge_merged, group.by = "old_celltype_completed", label = T, repel = T, split.by = "species", ncol = 3)

lunge_merged@active.ident = lunge_merged$SCT_snn_res.0.4

toolboxH::plotSankey(lunge_merged@meta.data[,c('SCT_snn_res.0.3',"old_celltype_completed", 'SCT_snn_res.0.4')], spalte4color = "old_celltype_completed")

## check difficult clusters ----

annolunge = lunge_merged@meta.data %>% as.data.table(keep.rownames = T)
annolunge[SCT_snn_res.0.3 %in% c(17,9),.N , .(celltypeReported,dataset,SCT_snn_res.0.3)][order(SCT_snn_res.0.3,N)]


# speichern-----
Seurat::GetAssay(lunge_merged)
lunge_merged@active.assay
lunge_merged@active.assay = "SCT"
lunge_merged@assays
lunge_merged@assays$RNA
Seurat::GetAssay(lunge_merged,  assay = "SCT")
saveRDS(lunge_merged, file = here("results/s603_2_seurat_lunge_merged_21-12-11.RDS"))
#


sct_features = rownames(GetAssayData(lunge_merged[["SCT"]], slot = "counts"))
str(sct_features)
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

saveRDS(lunge_mergedRNA, file = here("results/s603_2_seurat_NOsct_BBrowser_lunge_merged_21-12-11.RDS"))

finalizeSkript()




require(toolboxH)
require(here)
require(rmarkdown)

finalizeSkript()
#
# render(here("scripts/s603_2_harmonize_and_integrate_objects.R"),encoding="UTF-8")
