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



# # CREATE RENAMING TABLE for integration ----
gene_anno_pre = fread(here("data/mart_export.txt_orthologuesHumanCentered.gz"))  # downloaded from https://www.ensembl.org/biomart/martview/6dcbc1911b1c91a71d1d4055902e3c39  06/2021, Ensembl 104 GRCh38.p13

gene_anno_pre

setnames(gene_anno_pre, c("Gene name", "Gene stable ID", 'Gene description', "Gene stable ID version" ),  c("Human gene name", "Human gene stable ID",'Human gene description', "Human gene stable ID version"))

# Limit to genes defined in humans and existing in the scRNA-Seq data
gene_anno = gene_anno_pre[`Human gene name` != ""]

gene_anno[, name_hamster := `Golden Hamster gene name` ]
gene_anno[, name_mensch := `Human gene name`]
gene_anno[, name_mouse := `Mouse gene name`]


gene_anno4 = gene_anno[(name_hamster %in%   rownames(hamster@assays$RNA)  ) |
                             name_mensch  %in%   c(rownames(human_travaglini @assays$RNA),rownames(human_charite@assays$RNA)) |
                             name_mouse %in% rownames(mouse@assays$RNA)
]


# check length of genename as sanity control

gene_anno4$name_mensch %>% str_length() %>% table
gene_anno4$name_mouse %>% str_length() %>% table
gene_anno4$name_hamster %>%  str_length() %>% table


# reduce dataset to Human gene names, taking care of dupicates

# choose for duplicated genes higher expressed gene
# subset(hamster, features = unique(gene_anno4$name_hamster))

hamster@active.assay = "RNA"
hamster=NormalizeData(hamster )
av_express_hamster = AverageExpression(hamster)
hh(av_express_hamster$RNA)
dim(av_express_hamster$RNA)

human_travaglini@active.assay = "RNA"
human_travaglini=NormalizeData(human_travaglini)
av_express_human_travaglini  = AverageExpression(human_travaglini )
hh(av_express_human_travaglini$RNA)
dim(av_express_human_travaglini$RNA)

human_charite@active.assay = "RNA"
human_charite=NormalizeData(human_charite)
av_express_human_charite = AverageExpression(human_charite)
hh(av_express_human_charite$RNA)
dim(av_express_human_charite$RNA)

mouse@active.assay = "RNA"
mouse=NormalizeData(mouse)
av_express_mouse = AverageExpression(mouse)
hh(av_express_mouse$RNA)
dim(av_express_mouse$RNA)

av_express_mouse_dt = data.table(name_mouse = rownames(av_express_mouse$RNA),
                                 max_in_any_mouse = apply(av_express_mouse$RNA, 1, max))

av_express_mouse_dt

av_express_hamster_dt = data.table(name_hamster = rownames(av_express_hamster$RNA),
                                   max_in_any_hamster = apply(av_express_hamster$RNA, 1, max))

av_express_hamster_dt

av_express_human_travaglini_dt = data.table(name_mensch = rownames(av_express_human_travaglini$RNA),
                                             max_in_any_mensch = apply(av_express_human_travaglini$RNA, 1, max))

av_express_human_travaglini_dt


av_express_human_charite_dt = data.table(name_mensch = rownames(av_express_human_charite$RNA),
                                         max_in_any_mensch = apply(av_express_human_charite$RNA, 1, max))

av_express_human_charite_dt

stopifnot(nrow(av_express_human_travaglini_dt[allDuplicatedEntries(name_mensch)])==0)
stopifnot(nrow(av_express_human_charite_dt[allDuplicatedEntries(name_mensch)])==0)


# Define a single name for each species, preferring higher expressed genes when duplicates exists for hamster
av_express_hamster_dt2 = merge(av_express_hamster_dt, unique(gene_anno4[,.(name_mensch, name_hamster)]), by = 'name_hamster', all.x = T, sort = F)

av_express_hamster_dt2[name_mensch ==""]  %>% assertr::verify(nrow(.)==0)

av_express_hamster_dt2[,max_in_any_human_charite := av_express_human_charite_dt[match_hk(av_express_hamster_dt2$name_mensch, av_express_human_charite_dt$name_mensch, makeunique = T, importcol = av_express_human_charite_dt$max_in_any_mensch),max_in_any_mensch]]


av_express_hamster_dt2[, name_harmon := ifelse(is.na(name_mensch), name_hamster, name_mensch)]



showNA(av_express_hamster_dt2[, -c('name_mensch', 'max_in_any_human_charite')], showAllNoNA = F)  %>% assertr::verify(sum(.$NAs)==0)

setorder(av_express_hamster_dt2, -max_in_any_hamster)
av_express_hamster_dt2[is.na(name_mensch)==F][allDuplicatedEntries(name_mensch)]

av_express_hamster_dt2[is.na(name_mensch)==F, hamster_duplicate2remove := ifelse(duplicated(name_mensch), T, F)]
av_express_hamster_dt2[is.na(name_mensch)==T, hamster_duplicate2remove :=F]


setorder(av_express_hamster_dt2, -max_in_any_human_charite, na.last = T)
av_express_hamster_dt2[is.na(name_mensch)==F][allDuplicatedEntries(name_hamster)]
av_express_hamster_dt2[is.na(name_mensch)==F, mensch_duplicate2remove := ifelse(duplicated(name_hamster), T, F)]
av_express_hamster_dt2[is.na(mensch_duplicate2remove)==T, mensch_duplicate2remove :=F]

av_express_hamster_dt2[,.N, .(hamster_duplicate2remove,mensch_duplicate2remove)]


# Define a single name for each species, preferring higher expressed genes when duplicates exists for mouse
#
av_express_mouse_dt2 = merge(av_express_mouse_dt, unique(gene_anno4[,.(name_mensch, name_mouse)]), by = 'name_mouse', all.x = T, sort = F)

av_express_mouse_dt2[,max_in_any_human_charite := av_express_human_charite_dt[match_hk(av_express_mouse_dt2$name_mensch, av_express_human_charite_dt$name_mensch, makeunique = T, importcol = av_express_human_charite_dt$max_in_any_mensch),max_in_any_mensch]]


av_express_mouse_dt2[, name_harmon := ifelse(is.na(name_mensch), name_mouse, name_mensch)]

showNA(av_express_mouse_dt2[, -c('name_mensch', 'max_in_any_human_charite')], showAllNoNA = F) %>% assertr::verify(sum(.$NAs)==0)

setorder(av_express_mouse_dt2, -max_in_any_mouse)
av_express_mouse_dt2[is.na(name_mensch)==F][allDuplicatedEntries(name_mensch)]
av_express_mouse_dt2[is.na(name_mensch)==F, mouse_duplicate2remove := ifelse(duplicated(name_mensch), T, F)]
av_express_mouse_dt2[is.na(name_mensch)==T, mouse_duplicate2remove :=F]


setorder(av_express_mouse_dt2, -max_in_any_human_charite, na.last = T)
av_express_mouse_dt2[is.na(name_mensch)==F][allDuplicatedEntries(name_mouse)]
av_express_mouse_dt2[is.na(name_mensch)==F, mensch_duplicate2remove := ifelse(duplicated(name_mouse), T, F)]
av_express_mouse_dt2[is.na(name_mensch)==T, mensch_duplicate2remove :=F]

av_express_human_travaglini_dt[allDuplicatedEntries(name_mensch)] %>% assertr::verify(nrow(.)==0)
av_express_human_charite_dt[allDuplicatedEntries(name_mensch)] %>% assertr::verify(nrow(.)==0)

av_express_mouse_dt2[,.N, .(mouse_duplicate2remove,mensch_duplicate2remove)]


qlist_unique1 = venn4(av_express_human_charite_dt$name_mensch,
                      av_express_human_travaglini_dt$name_mensch,
                      av_express_hamster_dt2[hamster_duplicate2remove ==F & mensch_duplicate2remove==F, name_harmon],
                      av_express_mouse_dt2[mouse_duplicate2remove ==F & mensch_duplicate2remove==F, name_harmon], mylabels = c("human_charite","human_travaglini ", "hamster", "mouse"), mytitle = "Gene overlap accros species")

# # EXTEND ANNOTATION TABLE TO DOCUMENT RENAMING for later reference, including genes that have no orthologue in humans----

# replace "" with NA to assign missings in ensembl ata
for(i in names(gene_anno)) {
  gene_anno[,(i) := ifelse(get(i) =="", NA, get(i))]
}
gene_anno[,origin := "ensmart orthologues"]

# do this for hamster
qlist12_hamster = venn2(gene_anno$name_hamster, rownames(hamster))
str(qlist12_hamster)
no_orthologues_hamster_ensmaug = grep("^ENSMAUG", qlist12_hamster$q3, value=T)
no_orthologues_hamster_name = grep("^ENSMAUG", qlist12_hamster$q3,invert = T, value=T)

gene_anno2 = rbind(gene_anno,
                   data.table(`Golden Hamster gene name` = no_orthologues_hamster_name,
                              name_hamster = no_orthologues_hamster_name,
                              origin = "in hamster, no orthologues"),
                   data.table(`Golden Hamster gene stable ID` = no_orthologues_hamster_ensmaug,
                              name_hamster = no_orthologues_hamster_ensmaug,
                              origin = "in hamster, no orthologues"),
                   fill= T)

qlist12b_hamster = venn3(gene_anno2$`Golden Hamster gene name`, gene_anno2$`Golden Hamster gene stable ID`, rownames(hamster@assays$RNA))
stopifnot(length(qlist12b_hamster$q7)==0)

# mouse
qlist12_mouse = venn2(gene_anno2$`Mouse gene name`,  rownames(mouse@assays$RNA))
str(qlist12_mouse)
no_orthologues_mouse_name = qlist12_mouse$q3

gene_anno3 = rbind(gene_anno2, data.table(`Mouse gene name` = no_orthologues_mouse_name,
                                          name_mouse = no_orthologues_mouse_name,
                                          origin = "in mouse, no orthologues"), fill= T)


# do this for human travaglinie
qlist12_human_travaglini  = venn2(gene_anno3$`Human gene name`,  rownames(human_travaglini @assays$RNA))
str(qlist12_human_travaglini )
no_orthologues_human_travaglini_name = qlist12_human_travaglini$q3

gene_anno4 = rbind(gene_anno3, data.table(`Human gene name` = no_orthologues_human_travaglini_name,
                                          name_mensch = no_orthologues_human_travaglini_name,
                                          origin = "in mensch, no orthologues"), fill= T)
gene_anno4


# do this for human charite
qlist12_human_charite = venn2(gene_anno4$`Human gene name`,  rownames(human_charite@assays$RNA))
str(qlist12_human_charite)
no_orthologues_human_charite_name = qlist12_human_charite$q3

gene_anno5 = rbind(gene_anno4, data.table(`Human gene name` = no_orthologues_human_charite_name,
                                          name_mensch = no_orthologues_human_charite_name,
                                          origin = "in mensch, no orthologues"), fill= T)
gene_anno5
qlist12b_human_charite = venn2(gene_anno5$`Human gene name`,  rownames(human_charite@assays$RNA))

# replace again "" with NA to assign missings in ensembl ata
for(i in names(gene_anno5)) {
  gene_anno5[,(i) := ifelse(get(i) =="", NA, get(i))]
}


# document duplicates and expression level as filter criterium
gene_anno5[, max_in_any_human_travaglini  :=av_express_human_travaglini_dt[match_hk(gene_anno5$name_mensch, av_express_human_travaglini_dt$name_mensch), max_in_any_mensch]]

gene_anno5[, max_in_any_human_charite :=av_express_human_charite_dt[match_hk(gene_anno5$name_mensch, av_express_human_charite_dt$name_mensch), max_in_any_mensch]]

gene_anno5[, mouse_duplicate2remove := name_mouse %in% av_express_mouse_dt2[mouse_duplicate2remove ==T , name_mouse]]

gene_anno5[, mouse_mensch_duplicate2remove := name_mensch %in% av_express_mouse_dt2[mensch_duplicate2remove ==T , name_mensch]]

gene_anno5[, max_in_any_mouse :=av_express_mouse_dt2[match_hk(gene_anno5$name_mouse, av_express_mouse_dt2$name_mouse, makeunique = T, importcol = av_express_mouse_dt2$max_in_any_mouse), max_in_any_mouse]]

gene_anno5[, hamster_duplicate2remove := name_hamster %in% av_express_hamster_dt2[hamster_duplicate2remove ==T , name_hamster]]

gene_anno5[, hamster_mensch_duplicate2remove := name_mensch %in% av_express_hamster_dt2[mensch_duplicate2remove ==T , name_mensch]]

gene_anno5[, max_in_any_hamster :=av_express_hamster_dt2[match_hk(gene_anno5$name_hamster, av_express_hamster_dt2$name_hamster, makeunique = T, importcol = av_express_hamster_dt2$max_in_any_hamster), max_in_any_hamster]]


# ## define stress genes ----
stressgene_names =c("Fosb", "Fos", "Jun", "Junb", "Jund", "Atf3", "Egr1", "Hspa1a",
"Hspa1b", "Hsp90ab1", "Hspa8", "Hspb1", "Ier3", "Ier2", "Btg1",
"Btg2", "Dusp1")
gene_anno5[,stressgene := `Mouse gene name` %in% stressgene_names]



# ## define a consensus name, even if no orthologue existed----
gene_anno5[,gene_used_in_consensus := ifelse(is.na(name_mensch)==F, name_mensch, ifelse(is.na(name_mouse)==F, name_mouse, ifelse(is.na(name_hamster)==F, name_hamster, ifelse(is.na(`Human gene name`)==F, `Human gene name`, `Human gene stable ID`))))]


# # RENAME SEURAT OBJECTS using consensus name----
source(here("scripts/Seurat.update.gene.symbols.HGNC.R")) # from https://github.com/vertesy/Seurat.utils/

# ## rename hamster and check again for duplicates----
renaming_hamster = gene_anno5[hamster_duplicate2remove ==F & hamster_mensch_duplicate2remove==F & name_hamster %in% rownames(hamster), .(name_hamster, gene_used_in_consensus)] %>% unique
renaming_hamster[allDuplicatedEntries(name_hamster)] %>% assertr::verify(nrow(.)==0)


setorder(gene_anno5,gene_used_in_consensus, -max_in_any_human_charite, -max_in_any_human_travaglini ,-max_in_any_hamster , na.last = T)

renaming_hamster = gene_anno5[hamster_duplicate2remove ==F & hamster_mensch_duplicate2remove==F & name_hamster %in% rownames(hamster), .(name_hamster, gene_used_in_consensus)] %>% unique
renaming_hamster[allDuplicatedEntries(name_hamster)] %>% assertr::verify(nrow(.)==0)


renaming_hamster[,hamster_duplicate2remove := duplicated(gene_used_in_consensus) ]
renaming_hamster[allDuplicatedEntries(gene_used_in_consensus)]
renaming_hamster[hamster_duplicate2remove==T]
gene_anno5[paste(name_hamster, gene_used_in_consensus)==renaming_hamster[hamster_duplicate2remove==T, paste(name_hamster, gene_used_in_consensus)],hamster_duplicate2remove:= T]

renaming_hamster = renaming_hamster[hamster_duplicate2remove==F]
renaming_hamster[allDuplicatedEntries(gene_used_in_consensus)] %>% assertr::verify(nrow(.)==0)

showNA(renaming_hamster)%>% assertr::verify(sum(.$NAs)==0)
renaming_hamster[gene_used_in_consensus=="NA"]

qlisthamster =venn3(rownames(hamster), renaming_hamster$name_hamster, gene_anno5$name_hamster)
str(qlisthamster)
stopifnot(all(renaming_hamster$name_hamster %in% rownames(hamster)))
hamster3=hamster[renaming_hamster$name_hamster,]
hamster3
grep("^NA$", rownames(hamster3))%>% assertr::verify(length(.)==0)

stopifnot(identical(dimnames(hamster3@assays$RNA)[[1]], renaming_hamster$name_hamster))

hamster3 = RenameGenesSeurat(obj = hamster3, newnames = renaming_hamster$gene_used_in_consensus)
grep("^NA$", rownames(hamster3))%>% assertr::verify(length(.)==0)

rownames(hamster3)[duplicated(rownames(hamster3))]%>% assertr::verify(length(.)==0)


# ## rename mouse and check again for duplicates----

setorder(gene_anno5,gene_used_in_consensus, -max_in_any_human_charite, -max_in_any_human_travaglini ,-max_in_any_mouse , na.last = T)

renaming_mouse = gene_anno5[mouse_duplicate2remove ==F & mouse_mensch_duplicate2remove==F & name_mouse %in% rownames(mouse), .(name_mouse, gene_used_in_consensus)] %>% unique
renaming_mouse[allDuplicatedEntries(name_mouse)] %>% assertr::verify(nrow(.)==0)

renaming_mouse[,mouse_duplicate2remove := duplicated(gene_used_in_consensus) ]
renaming_mouse[allDuplicatedEntries(gene_used_in_consensus)]
renaming_mouse[mouse_duplicate2remove==T]
gene_anno5[paste(name_mouse, gene_used_in_consensus)==renaming_mouse[mouse_duplicate2remove==T, paste(name_mouse, gene_used_in_consensus)],mouse_duplicate2remove:= T]

renaming_mouse = renaming_mouse[mouse_duplicate2remove==F]


renaming_mouse[allDuplicatedEntries(gene_used_in_consensus)] %>% assertr::verify(nrow(.)==0)

showNA(renaming_mouse)%>% assertr::verify(sum(.$NAs)==0)
renaming_mouse[gene_used_in_consensus=="NA"]

qlistmouse =venn3(rownames(mouse), renaming_mouse$name_mouse, gene_anno5$name_mouse)
stopifnot(all(renaming_mouse$name_mouse %in% rownames(mouse)))
mouse3=mouse[renaming_mouse$name_mouse,]
mouse3
grep("^NA$", rownames(mouse3))%>% assertr::verify(length(.)==0)

stopifnot(identical(dimnames(mouse3@assays$RNA)[[1]], renaming_mouse$name_mouse))

mouse3 = RenameGenesSeurat(obj = mouse3, newnames = renaming_mouse$gene_used_in_consensus)
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


# # SAVE ----

fwrite(gene_anno5, file =here("results/s1_1_gene_annotation_orthologues_harmonized4integration.txt"), sep = "\t")

saveRDS(lung_list, here("results/s1_1_list_human_mouse_hamster_harmonized4integration.RDS"))


# # STATS
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

# # FINALIZE----

finalizeSkript()


