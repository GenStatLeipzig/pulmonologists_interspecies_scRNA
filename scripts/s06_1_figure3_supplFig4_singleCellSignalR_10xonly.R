rm(list = ls())

library(SingleCellSignalR)
library(Seurat)
library(generics)
library(here)
library(dplyr)
require(toolboxH)
require(stringr)


# # LOAD seurat objekt laden----
seurat_fn = here("results/s603_2_seurat_lunge_merged_21-12-11.RDS")
stopifnot(file.exists(seurat_fn))
seurat_fn
seurat = readRDS(seurat_fn)

uniqueN(rownames(seurat))

# apply current cell type names ----
seurat$celltype <- paste("C", seurat$SCT_snn_res.0.3)
seurat@active.ident = factor(seurat$celltype)

seurat <- RenameIdents (seurat, "C 0" = "T Cells", 'C 1' = "Endothelial", 'C 2' = "AM", 'C 3' = "Macrophages/ Monocytes", 'C 4' = "NK Cells", 'C 5' = "AT2", 'C 6' = "Endothelial", 'C 7' = "Fibroblasts", 'C 8' = "B Cells", 'C 9' = "Perivascular", 'C 10' = "AT1", 'C 11' = "Macrophages/ Monocytes", 'C 12' = "Mast Cells", 'C 13' = "Ciliated Cells", 'C 14' = "Club Cells", 'C 15' = "Dendritic Cells", 'C 16' = "Macrophages/ Monocytes", 'C 17' = "ly Endothelial", 'C 18' = "Proliferating NK/T", "C 19" = "Proliferating AM", "C 20" = "Neutrophils", "C 21" = "Endothelial", "C 22" = "Fibroblasts")


seurat$celltype_pre = Idents(seurat)
seurat$celltype = str_replace_all(seurat$celltype_pre, " |/", "_")

seurat = SetIdent(seurat, value = "celltype")

# LOAD Ligand - REceptor pairs from Raredon 2019 Supplement ----
raredon = read_excel2(here("data/raredon_2019_pmid31840053_Single-cell_connectomic_analysis_of_adult_mammalia_suppl1.xlsx"))
raredon[, receptor := str_split(pair, "-") %>% sapply(., "[", 2)]
raredon[, ligand := str_split(pair, "-") %>% sapply(., "[", 1)]
showNA(raredon)

# CHeck overlap with data.base SingleCellSignalR ----
LRdb2 = as.data.table(LRdb, keep.rownames = T)
LRdb2[, pair:= paste0(ligand, "-", receptor)]

qlist234 = venn2(raredon$pair, LRdb2$pair)
raredon[pair %in% qlist234$q2]


qlist1 = venn3(raredon$receptor, raredon$ligand, rownames(seurat)) # all included

candidate_genes = rbind(raredon[, .(typ = mode, hgnc = ligand)], 
                        raredon[, .(typ = mode, hgnc = receptor)]
)

# # REduce Seurat Object to Raredon celltypes and 10x ----
table(Idents(seurat))
raredon_celltypes = c("Endothelial",
                      "AM",
                      "AT2",
                      "Perivascular",
                      "AT1",
                      "Proliferating_AM",
                      "Fibroblasts")
qlist33 =venn2(Idents(seurat), raredon_celltypes)
stopifnot(length(qlist33$q3)==0)

table(seurat$dataset)
seurat4 = seurat[,Idents(seurat) %in% raredon_celltypes &
                   seurat$dataset != 'humanSS2' ]
table(seurat4$dataset)
seurat4

DimPlot(seurat4,label =  T ) + NoLegend()



# # >Ligand/Receptor analysis using SingleCellSignalR berechnen-----

calculateInteractions = function(seuratojekt, 
                                 LRScore_cutoff=0.5, # for allternate scenario
                                 slot2use = "SCT",
                                 doscaleGenes = F) {
  # seuratojekt = seurat4
  # LRScore_cutoff=0.5
  # slot2use = "SCT"
  # doscaleGenes=F
  
  cluster_numeric = as.numeric(Idents(seuratojekt))
  # table(cluster_numeric)
  
  data_RNA = data.frame(seuratojekt[[slot2use]]@data) # SCT or RNA 33944 x 85098 sparse Matrix of class "dgCMatrix" wenn all
  hh(data_RNA)
  all.genes = rownames(seuratojekt)
  if(doscaleGenes) {
    warning("Scaling experimental, not meaningful for standard - analyses")
    data_RNA = apply(data_RNA, 2, scale) 
    rownames(data_RNA) = all.genes
    data_RNA[is.nan(data_RNA)] = 0
  }
  # str(all.genes)
  
  
  # cell-cell-interakt-alle----
  cluster_name2num = data.table(name = unique(Idents(seuratojekt)))
  cluster_name2num[,num := as.numeric(name)]
  cluster_name2num
  
  message("Calculating Interaction without cutoff (LRscore = 0)")
  interactions_10x = c()
  
  hh(data_RNA)
  interactions_10x$singleCellSignalR_all = cell_signaling(data=data_RNA,genes=all.genes,cluster=cluster_numeric,c.names = cluster_name2num[order(num), paste0(as.character(name), "\n")], write = F, s.score = 0)
  # interactions_10x$singleCellSignalR_all
  
  message("Calculating Interaction with cutoff (LRscore = ", LRScore_cutoff,")")
  interactions_10x$singleCellSignalR_cutoff_0.5 = cell_signaling(data=data_RNA,genes=all.genes,cluster=cluster_numeric,c.names = cluster_name2num[order(num), paste0(as.character(name), "\n")], write = F, s.score = LRScore_cutoff)
  # interactions_10x$singleCellSignalR_cutoff_0.5
  #
  #
  # signal_peter2_05auto = cell_signaling(data=data_RNA,genes=all.genes,cluster=cluster_numeric,c.names = cluster_name2num[order(num), paste0(as.character(name), "\n")], write = F, s.score = 0.5,int.type = c( "autocrine"))
  # signal_peter2_05auto
  # visualize_interactions(signal_peter2_05auto)
  
  
  
  signal_peter_detail = lapply(interactions_10x$singleCellSignalR_all, function(ll) {
    # ll = signal_peter2$`cluster 1-cluster 2`
    #ll= interactions_10x$singleCellSignalR_all[1]
    ll = data.table(ll)
    ll[,cluster_a:= names(ll)[1]]
    ll[,cluster_b:= names(ll)[2]]
    names(ll)[1] = "gene_a"
    names(ll)[2] = "gene_b"
    ll
  }) %>% rbindlist
  
  signal_peter_detail
  
  chec1 = venn2(signal_peter_detail$gene_a, signal_peter_detail$gene_b, plotte = F)
  chec1$q1
  # signal_peter_detail[gene_a %in% chec1$q1 | gene_b %in% chec1$q1]
  # signal_peter_detail[,.N,`interaction type` ]
  #
  signal_peter_detail[,cluster_a_name := str_replace(cluster_a, "\\n", "")]
  signal_peter_detail[,cluster_b_name := str_replace(cluster_b, "\\n", "")]
  
  res = c()
  res$singleCellSignalR_all = interactions_10x$singleCellSignalR_all
  
  res$singleCellSignalR_cutoff_0.5 = interactions_10x$singleCellSignalR_cutoff_0.5
  res$singleCellSignalR_all_table = signal_peter_detail
  res
}


interactions_10x = calculateInteractions(seuratojekt = seurat4) # warnings "No such file as table_dge_T Cells .txt in the cluster-analysis folder" means no differential expression done. This is an optional step to characterize interactions in diff. expr. genes, only
visualize_interactions(interactions_10x$singleCellSignalR_all)
visualize_interactions(interactions_10x$singleCellSignalR_cutoff_0.5)


savePlotInMultiFormats(here("results/s606_9_signalCellSignalR_overview_limit0.5_ACX"), weite =  6,hoehe =  6,plotecode = 'visualize_interactions(interactions_10x$singleCellSignalR_cutoff_0.5)')




# # Now calculate species-specific ---- 

interactions_human = calculateInteractions(seuratojekt = seurat4[,grepl("human", seurat4$dataset)]) 


savePlotInMultiFormats(here("results/s606_9_signalCellSignalR_humanALL_limit0.5_ACX"), weite =  6,hoehe =  6,plotecode = 'visualize_interactions(interactions_human$singleCellSignalR_cutoff_0.5)')

visualize_interactions(interactions_human$singleCellSignalR_cutoff_0.5)

interactions_monkey = calculateInteractions(seuratojekt = seurat4[,grepl("monkey", seurat4$dataset)]) 
visualize_interactions(interactions_monkey$singleCellSignalR_cutoff_0.5)


names(interactions_human$singleCellSignalR_cutoff_0.5)
names(interactions_monkey$singleCellSignalR_cutoff_0.5)

interactions_pig = calculateInteractions(seuratojekt = seurat4[,grepl("pig", seurat4$dataset)]) 
visualize_interactions(interactions_pig$singleCellSignalR_cutoff_0.5)

interactions_hamster = calculateInteractions(seuratojekt = seurat4[,grepl("hamster", seurat4$dataset)]) 
visualize_interactions(interactions_hamster$singleCellSignalR_cutoff_0.5)

interactions_mouse = calculateInteractions(seuratojekt = seurat4[,grepl("mouse", seurat4$dataset)]) 
visualize_interactions(interactions_mouse$singleCellSignalR_cutoff_0.5)

interactions_rat = calculateInteractions(seuratojekt = seurat4[,grepl("rat", seurat4$dataset)]) 
visualize_interactions(interactions_rat$singleCellSignalR_cutoff_0.5)

pdf(here("results/s606_9_human_chordplot.pdf"), 6, 6)
visualize_interactions(interactions_human$singleCellSignalR_cutoff_0.5)
title("HUMAN", cex.main=1.5)
dev.off()

pdf(here("results/s606_9_hamster_chordplot.pdf"), 6, 6)
visualize_interactions(interactions_hamster$singleCellSignalR_cutoff_0.5)
title("HAMSTER", cex.main=1.5)
dev.off()

pdf(here("results/s606_9_mouse_chordplot.pdf"), 6, 6)
visualize_interactions(interactions_mouse$singleCellSignalR_cutoff_0.5)
title("MOUSE", cex.main=1.5)
dev.off()

pdf(here("results/s606_9_pig_chordplot.pdf"), 6, 6)
visualize_interactions(interactions_pig$singleCellSignalR_cutoff_0.5)
title("PIG", cex.main=1.5)
dev.off()

pdf(here("results/s606_9_monkey_chordplot.pdf"), 6, 6)

visualize_interactions(interactions_monkey$singleCellSignalR_cutoff_0.5)
title("MONKEY", cex.main=1.5)

dev.off()

pdf(here("results/s606_9_rat_chordplot.pdf"), 6, 6)
visualize_interactions(interactions_rat$singleCellSignalR_cutoff_0.5)
title("RAT", cex.main=1.5)
dev.off()

pdf(here("results/s606_9_all_species_chordplot.pdf"), 6, 6)
visualize_interactions(interactions_human$singleCellSignalR_cutoff_0.5)
title("HUMAN", cex.main=1.5)
visualize_interactions(interactions_hamster$singleCellSignalR_cutoff_0.5)
title("HAMSTER", cex.main=1.5)
visualize_interactions(interactions_mouse$singleCellSignalR_cutoff_0.5)
title("MOUSE", cex.main=1.5)
visualize_interactions(interactions_pig$singleCellSignalR_cutoff_0.5)
title("PIG", cex.main=1.5)
visualize_interactions(interactions_monkey$singleCellSignalR_cutoff_0.5)
title("MONKEY", cex.main=1.5)
visualize_interactions(interactions_rat$singleCellSignalR_cutoff_0.5)
title("RAT", cex.main=1.5)
dev.off()
# 
makeFocusedTable <- function(candidate_genes, interactions_10x, subexperiment="", LAbel_interCategs_as = "Other", receptor_ligand = "und", no_between_group_interaction=F , match_or_merge = "match", LRScore_cutoff=0.5, group_plot_by = "ligand", custom_x_order = NULL, nrow_legend=2) {
  
  # match_or_merge = "merge"
  showplots = F
  # ueberlapp
  qlist5 = venn3(candidate_genes$hgnc, interactions_10x$singleCellSignalR_all_table$gene_a, interactions_10x$singleCellSignalR_all_table$gene_b, plotte = showplots)
  
  candidate_genes
  
  
  if(receptor_ligand=="oder")   signal_peter_detail_fokus = interactions_10x$singleCellSignalR_all_table[gene_a %in% na.omit(candidate_genes$hgnc) |
                                                                                                           gene_b %in% na.omit(candidate_genes$hgnc)] else if(receptor_ligand=="und")  signal_peter_detail_fokus = interactions_10x$singleCellSignalR_all_table[gene_a %in% na.omit(candidate_genes$hgnc) & gene_b %in% na.omit(candidate_genes$hgnc)] else stop("parameter 'receptor_ligand' muss `oder` bzw. `und` sein")
  
  signal_peter_detail_fokus[, clustpair := paste0(cluster_a_name, " -> ", cluster_b_name)]
  signal_peter_detail_fokus[, genepair := paste0(gene_a, " / ", gene_b)]
  
  if(match_or_merge == "match") {
    signal_peter_detail_fokus[,gruppe_a := candidate_genes[match_hk(signal_peter_detail_fokus$gene_a, candidate_genes$hgnc, makeunique = T, importcol = candidate_genes$typ), typ]]
    
    signal_peter_detail_fokus[,gruppe_b := candidate_genes[match_hk(signal_peter_detail_fokus$gene_b, candidate_genes$hgnc, makeunique = T, importcol = candidate_genes$typ), typ]]
    
  }
  
  if(match_or_merge == "merge") {
    
    signal_peter_detail_fokus = merge(signal_peter_detail_fokus %>% unique(), candidate_genes[,.(hgnc, gruppe_a=typ)] %>% unique(), by.x = "gene_a", by.y = "hgnc", all.x = T, allow.cartesian=TRUE)
    
    signal_peter_detail_fokus = merge(signal_peter_detail_fokus %>% unique(), candidate_genes[,.(hgnc, gruppe_b=typ)] %>% unique(), by.x = "gene_b", by.y = "hgnc", all.x = T, allow.cartesian=TRUE)
    
    signal_peter_detail_fokus = unique(signal_peter_detail_fokus)
  } else stop('match_or_merge must be either "merge" or  "match"')
  
  if(no_between_group_interaction == T) signal_peter_detail_fokus = signal_peter_detail_fokus[gruppe_a ==gruppe_b]
  
  signal_peter_detail_fokus[,zeile := .I]
  signal_peter_detail_fokus[, gruppe := unique(na.omit(c(gruppe_a, gruppe_b))) %>% paste(., collapse = "+"), zeile] 
  
  if(length(LAbel_interCategs_as)>0) {
    signal_peter_detail_fokus[grep("\\+", gruppe), gruppe := LAbel_interCategs_as]
  }
  signal_peter_detail_fokus[,.N, gruppe]
  
  candidate_genes[,in_SingleCellSignalR_resultAll := hgnc %in% interactions_10x$singleCellSignalR_all_table[, c(gene_a, gene_b)]]
  
  candidate_genes[in_SingleCellSignalR_resultAll==T]
  
  
  require(ggplot2)
  require(scales)
  require(ggthemes)
  
  signal_peter_detail_fokus[gruppe =="Tumor Nekrose Faktor", gruppe := "Tumor\nNekrose\nFaktor"]
  
  
  ## add 0 if no value present
  signal_peter_detail_fokus[,max_LRScore := max(LRscore, na.rm = T),genepair  ]
  allcombis = expand.grid(clustpair = signal_peter_detail_fokus[max_LRScore>=LRScore_cutoff, unique(clustpair)],
                          genepair = signal_peter_detail_fokus[max_LRScore>=LRScore_cutoff, unique(genepair)]) %>% data.table()
  allcombis[,gruppe := signal_peter_detail_fokus[match_hk(allcombis$genepair, signal_peter_detail_fokus$genepair, makeunique = T, importcol = signal_peter_detail_fokus$gruppe), gruppe] ]
  
  allcombis[,cluster_a_name := str_split(clustpair, "->") %>% sapply(., "[",1) %>% str_trim]
  allcombis[,cluster_b_name := str_split(clustpair, "->") %>% sapply(., "[",2) %>% str_trim]
  
  allcombis2 = allcombis[paste(clustpair, genepair)%nin% signal_peter_detail_fokus[,paste(clustpair, genepair)]]
  allcombis2[,LRscore:=0]
  
  signal_peter_detail_fokus2 = rbind(signal_peter_detail_fokus, allcombis2, fill = T)
  signal_peter_detail_fokus2[,max_LRScore := max(LRscore, na.rm = T),genepair  ]
  signal_peter_detail_fokus2[,clustpair := factor(clustpair, levels = sort(unique(as.character(clustpair))))]
  
  setorder(signal_peter_detail_fokus2, max_LRScore)
  signal_peter_detail_fokus2[,genepair := factor(genepair, levels = unique(as.character(genepair)))]
  
  signal_peter_detail_fokus2[,max_LRScore_row := max(LRscore, na.rm = T),clustpair  ]
  setorder(signal_peter_detail_fokus2, -max_LRScore_row)
  signal_peter_detail_fokus2[,clustpair := factor(clustpair, levels = unique(as.character(clustpair)))]
  
  signal_peter_detail_fokus2$subexperiment = subexperiment
  
  if(is.null(custom_x_order)==F) {
    signal_peter_detail_fokus2[, clustpair := factor(clustpair, levels = sort(unique(clustpair))[custom_x_order])]
    message("using levels ", paste(paste(1:length(levels(signal_peter_detail_fokus2$clustpair)), levels(signal_peter_detail_fokus2$clustpair)), collapse= "\n"))
  }
  
  p_fokus = ggplot(signal_peter_detail_fokus2[max_LRScore>=LRScore_cutoff & max_LRScore_row>=LRScore_cutoff], aes(clustpair, y= genepair, fill= -LRscore, label = round(10*LRscore) )) +
    geom_tile(   colour = "white") + 
    facet_grid(gruppe~cluster_a_name, scales = "free", space = "free" ) +scale_fill_gradient2_tableau(labels = c(0.8,0.7, 0.6, 0.5,0.4, 0.3,0.1,0), breaks = -c(0.8,0.7, 0.6, 0.5,0.4, 0.3,0.1,0) , guide = "legend") + 
    labs(fill = "LR score" )+ 
    guides(fill=guide_legend(nrow=1))+
    # theme_light(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
          strip.text.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text.y = element_text(color = "black", face = "bold", angle  =0),
          strip.background.y = element_rect(fill = "white"),
          legend.position = "top"
    ) + xlab("") # + coord_fixed(ratio=1)
  
  if(group_plot_by == "receptor") {
    
    p_fokus = p_fokus +  facet_grid(gruppe~cluster_b_name, scales = "free", space = "free" )
    
    
  }
  
  
  
  plot(p_fokus)
  
  # tabelle nach auffuellen annotieren
  
  signal_peter_detail_fokus2[,gene_a := str_split(genepair, " / ") %>% sapply(., "[",1) %>% str_trim]
  signal_peter_detail_fokus2[,gene_b := str_split(genepair, " / ") %>% sapply(., "[",2) %>% str_trim]
  signal_peter_detail_fokus2[,gene_a_in_cand:= gene_a %in% candidate_genes$hgnc]
  signal_peter_detail_fokus2[,gene_b_in_cand:= gene_b %in% candidate_genes$hgnc]
  
  # print(signal_peter_detail_fokus2)
  
  res = c()
  res$plot = p_fokus
  res$table = signal_peter_detail_fokus2[max_LRScore>=LRScore_cutoff &  max_LRScore_row>=LRScore_cutoff, .(subexperiment, genepair, clustpair, LRscore, max_LRScore, gruppe,cluster_a_name,cluster_b_name, gene_a, gene_b, gene_a_in_cand, gene_b_in_cand)]
  
  
  res$table_all = signal_peter_detail_fokus2
  
  res
}




focusobjekt_10x = makeFocusedTable(candidate_genes , interactions_10x, subexperiment = "10x", match_or_merge = "merge", no_between_group_interaction = T)
focusobjekt_10x$plot

focusobjekt_human = makeFocusedTable(candidate_genes , interactions_human, subexperiment = "human", match_or_merge = "merge", no_between_group_interaction = T)
focusobjekt_human$plot


focusobjekt_pig = makeFocusedTable(candidate_genes , interactions_pig, subexperiment = "pig", match_or_merge = "merge", no_between_group_interaction = T) 

focusobjekt_rat = makeFocusedTable(candidate_genes , interactions_rat, subexperiment = "rat", match_or_merge = "merge", no_between_group_interaction = T) 

focusobjekt_monkey = makeFocusedTable(candidate_genes , interactions_monkey, subexperiment = "monkey", match_or_merge = "merge", no_between_group_interaction = T) 

focusobjekt_hamster = makeFocusedTable(candidate_genes , interactions_hamster, subexperiment = "hamster", match_or_merge = "merge", no_between_group_interaction = T) 

focusobjekt_mouse = makeFocusedTable(candidate_genes , interactions_mouse, subexperiment = "mouse", match_or_merge = "merge", no_between_group_interaction = T) 

# conserved ----
# now plot genepairs following the classification of Raredon et al.

allspeciesFocusTab = rbind(focusobjekt_mouse$table_all %>% mutate(species = "mouse"),
                           focusobjekt_rat$table_all %>% mutate(species = "rat"),
                           focusobjekt_monkey$table_all %>% mutate(species = "monkey"),
                           focusobjekt_pig$table_all %>% mutate(species = "pig"),
                           focusobjekt_hamster$table_all %>% mutate(species = "hamster"),
                           focusobjekt_human$table_all %>% mutate(species = "human")
)
                           
                           
allspeciesFocusTab[,.N, .(gene_a_in_cand, gene_b_in_cand)]
stopifnot(b=nrow(allspeciesFocusTab[LRscore >=0.5 & gruppe_a != gruppe_b])==0)
focusobjekt_human_raredon= focusobjekt_human$table[(paste0(gene_a, "-",gene_b) %in% raredon$pair) & (cluster_a_name != cluster_b_name)] # this considers all receptor - ligand interaction based on the pair-wise classification of Raredon
# focusobjekt_human_raredon = copy(focusobjekt_human$table_all) # this would consider all receptor - ligand interaction based on the gene-wise classification of Raredon



qlist4 = venn2(focusobjekt_human_raredon[LRscore>=0.5 ,genepair], allspeciesFocusTab[LRscore>=0.5 & species != "human", genepair])
qlist5 = venn2(focusobjekt_human_raredon[LRscore>=0.5 ,clustpair], allspeciesFocusTab[LRscore>=0.5 & species != "human", clustpair])

focusobjekt_human_raredon_conserved = focusobjekt_human_raredon[genepair %in% qlist4$q1 & clustpair %in% qlist5$q1]

focusobjekt_human_raredon_conserved
focusobjekt_human_raredon_conserved[,max_LRScore := max(LRscore, na.rm = T),genepair  ]
focusobjekt_human_raredon_conserved[,clustpair := factor(clustpair, levels = sort(unique(as.character(clustpair))))]

setorder(focusobjekt_human_raredon_conserved, max_LRScore)
focusobjekt_human_raredon_conserved[,genepair := factor(genepair, levels = unique(as.character(genepair)))]

focusobjekt_human_raredon_conserved[,max_LRScore_row := max(LRscore, na.rm = T),clustpair  ]
setorder(focusobjekt_human_raredon_conserved, -max_LRScore_row)
focusobjekt_human_raredon_conserved[,clustpair := factor(clustpair, levels = unique(as.character(clustpair)))]


LRScore_cutoff = 0.5
p_fokus_human_raredon_conserved = ggplot(focusobjekt_human_raredon_conserved[max_LRScore>=LRScore_cutoff & max_LRScore_row>=LRScore_cutoff], aes(clustpair, y= genepair, fill= -LRscore, label = round(10*LRscore) )) +
  geom_tile(   colour = "white") + 
  facet_grid(gruppe~cluster_a_name, scales = "free", space = "free" ) +scale_fill_gradient2_tableau(labels = c(0.8,0.7, 0.6, 0.5,0.4, 0.3,0.1,0), breaks = -c(0.8,0.7, 0.6, 0.5,0.4, 0.3,0.1,0) , guide = "legend") + 
  labs(fill = "LR score" )+ 
  guides(fill=guide_legend(nrow=1))+
  # theme_light(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(color = "black", face = "bold", angle  =0),
        strip.background.y = element_rect(fill = "white"),
        legend.position = "top"
  ) + xlab("")

p_fokus_human_raredon_conserved

pdf(here("results/s606_9_conservedONLY_raredon_interactions.pdf"), 6,9.5)
p_fokus_human_raredon_conserved
dev.off()



# # Now plot supplementary table incl. uman interaction ----

focusobjekt_human_raredon
focusobjekt_human_raredon[,max_LRScore := max(LRscore, na.rm = T),genepair  ]
focusobjekt_human_raredon[,clustpair := factor(clustpair, levels = sort(unique(as.character(clustpair))))]

setorder(focusobjekt_human_raredon, max_LRScore)
focusobjekt_human_raredon[,genepair := factor(genepair, levels = unique(as.character(genepair)))]

focusobjekt_human_raredon[,max_LRScore_row := max(LRscore, na.rm = T),clustpair  ]
setorder(focusobjekt_human_raredon, -max_LRScore_row)
focusobjekt_human_raredon[,clustpair := factor(clustpair, levels = unique(as.character(clustpair)))]


LRScore_cutoff = 0.5
focusobjekt_human_raredon[gruppe =="Vasoactive" & LRscore>=0.5][order(genepair)]

focusobjekt_human_raredon_plot = focusobjekt_human_raredon[max_LRScore>=LRScore_cutoff & max_LRScore_row>=LRScore_cutoff]
focusobjekt_human_raredon_plot[gruppe =="Vasoactive" & LRscore>=0.5, ][order(genepair)]

p_focusobjekt_human_raredon = ggplot(focusobjekt_human_raredon_plot, aes(clustpair, y= genepair, fill= -LRscore, label = round(10*LRscore) )) +
  geom_tile(   colour = "white") + 
  facet_grid(gruppe~cluster_a_name, scales = "free", space = "free" ) +scale_fill_gradient2_tableau(labels = c(0.8,0.7, 0.6, 0.5,0.4, 0.3,0.1,0), breaks = -c(0.8,0.7, 0.6, 0.5,0.4, 0.3,0.1,0) , guide = "legend") + 
  labs(fill = "LR score" )+ 
  guides(fill=guide_legend(nrow=1))+
  # theme_light(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(color = "black", face = "bold", angle  =0),
        strip.background.y = element_rect(fill = "white"),
        legend.position = "top"
  ) + xlab("")

p_focusobjekt_human_raredon

# add stars for conserved pairs in any of eachs interactions ----


qlist8 = venn2(focusobjekt_human_raredon_conserved[LRscore>0.5, genepair], focusobjekt_human_raredon[LRscore>0.5, genepair])
qlist9 = venn2(focusobjekt_human_raredon_conserved[LRscore>0.5, clustpair], focusobjekt_human_raredon[LRscore>0.5, clustpair])

focusobjekt_human_raredon_plot[,genepair2 :=  ifelse(genepair %in% qlist8$q1, paste("*", genepair), genepair %>% as.character)]
focusobjekt_human_raredon_plot[,clustpair2 :=  ifelse(clustpair %in% qlist9$q1, paste("*", clustpair),clustpair %>% as.character)]

p_focusobjekt_human_raredon_conservedinfo = ggplot(focusobjekt_human_raredon_plot, aes(clustpair2, y= genepair2, fill= -LRscore, label = round(10*LRscore) )) +
  geom_tile(   colour = "white") + 
  facet_grid(gruppe~cluster_a_name, scales = "free", space = "free" ) +scale_fill_gradient2_tableau(labels = c(0.8,0.7, 0.6, 0.5,0.4, 0.3,0.1,0), breaks = -c(0.8,0.7, 0.6, 0.5,0.4, 0.3,0.1,0) , guide = "legend") + 
  labs(fill = "LR score" )+ 
  guides(fill=guide_legend(nrow=1))+
  # theme_light(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(color = "black", face = "bold", angle  =0),
        strip.background.y = element_rect(fill = "white"),
        legend.position = "top"
  ) + xlab("")

p_focusobjekt_human_raredon_conservedinfo

pdf(here("results/s606_9_conserved_raredon_interactions.pdf"), 6.6,12)
p_focusobjekt_human_raredon_conservedinfo
dev.off()

finalizeSkript()
