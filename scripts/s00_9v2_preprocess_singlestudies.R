rm(list = ls())
require(toolboxH)
require(here)
require(Seurat)
require(ggplot2)
require(ggthemes)
require(dplyr)
require(tidyr)
require(patchwork)


# Load gene annotation accross species ----
orthologues = fread(here("data/martquery_1119125514_551.txt.gz")) # "Ensembl release 104 - May 2021 Human genes (GRCh38.p13), https://www.ensembl.org/biomart/martview/bf27d7b2470fe6d55826530cdb586229

setnames(orthologues, c("Gene name", "Gene stable ID", 'Gene description' ),  c("Human gene name", "Human gene stable ID",'Human gene description'))

param_QC_filter_mito_max 	=30 # ABX 30, human 10
param_QC_filter_RNAfeatures_min 	=300
param_QC_filter_RNAfeatures_max 	= 6000
param_QC_filter_RNAcounts_min 	=1000
param_QC_filter_RNAcounts_max	=35000

# replace "" with NA in Gene names----
showNA(orthologues, showAllNoNA = F)
for(i in c("Human gene name",  "Mouse gene name",'Human gene stable ID', 'Mouse gene stable ID')) {
  orthologues[get(i)=="",(i):=NA]
}
showNA(orthologues, showAllNoNA = F)


# # >> FUNCTIONS---
plotQC = function(seuratobject) {
  qccols = c('run10x', 'nCount_RNA', 'nFeature_RNA',  'pct.mito',
             'pct.ribo')
  plotdata = seuratobject@meta.data %>% dplyr::select(qccols[qccols %in% names(seuratobject@meta.data)])
  p1 = ggplot( plotdata %>%
                 tidyr::pivot_longer(-run10x, names_to = 'metric', values_to = 'value'),
               aes(x=run10x,y=value,fill=run10x)) +
    geom_boxplot(outlier.size=.5) +
    facet_wrap(~metric,ncol=4,scales='free_y') +
    theme(axis.text.x=element_blank())

  if('pct.mito' %in% names(seuratobject@meta.data)) {
    p2 = ggplot(plotdata, aes(x=nCount_RNA ,y=pct.mito)) +
      geom_point(alpha = 0.3) +
      facet_wrap(~run10x,ncol=4,scales='free_y')

    p3 = ggplot(plotdata, aes(x=nFeature_RNA ,y=pct.mito)) +
      geom_point(alpha = 0.3) +
      facet_wrap(~run10x,ncol=4,scales='free_y')

    p4 = ggplot(plotdata, aes(x=nCount_RNA ,y=nFeature_RNA)) +
      geom_point(alpha = 0.3) +
      facet_wrap(~run10x,ncol=4,scales='free_y')

    plot(p1+p2+p3+p4 + patchwork::plot_layout(ncol = 1))
  } else {


    p4 = ggplot(plotdata, aes(x=nCount_RNA ,y=nFeature_RNA)) +
      geom_point(alpha = 0.3) +
      facet_wrap(~run10x,ncol=4,scales='free_y')

    plot(p1+p4 + patchwork::plot_layout(ncol = 1))
  }

}

doBAsicFiltering = function(seuratobject,
                            QC_filter_mito_max=param_QC_filter_mito_max,
                            QC_filter_RNAfeatures_max=param_QC_filter_RNAfeatures_max,
                            QC_filter_RNAfeatures_min=param_QC_filter_RNAfeatures_min,
                            QC_filter_RNAcounts_max=param_QC_filter_RNAcounts_max,
                            QC_filter_RNAcounts_min=param_QC_filter_RNAcounts_min) {
  # seuratobject= monkey
  seuratobject_name = deparse(substitute(seuratobject))
  dim_before = dim(seuratobject)
  if("pct.mito" %in% names(seuratobject@meta.data)) {
    seuratobject <- subset(seuratobject, (pct.mito < QC_filter_mito_max))
  }


  seuratobject <- subset(seuratobject,
                         (nFeature_RNA < QC_filter_RNAfeatures_max) &
                           (nFeature_RNA > QC_filter_RNAfeatures_min) &
                           (nCount_RNA < QC_filter_RNAcounts_max) &
                           (nCount_RNA > QC_filter_RNAcounts_min) )
  dim_after = dim(seuratobject)
  message(Sys.time(),"...Filtering ", seuratobject_name, " from\n", dim_before[1], " genes x ", dim_before[2], " cells to \n", dim_after[1], " genes x ", dim_after[2], " cells")
  message("(i.e. removing ", ((dim_before[2]-dim_after[2])/dim_after[2]) %>% proz()," of given cells)")
  seuratobject
}


correctAmbientRNA_inclEmptyCells = function(seuratobject, raw_h5_10x_fn, experiment_name,empty_drops_lower = 100, ribogenes=NULL, mitogenes=NULL ,
                                            orthologue_genenames = NULL, # if mitogenes or ribogenes are provided than  if also providing orthologue_genenames, percentage of ribogenes and mitogenes relates to orthologue_genenames also included in the seuratobject, only
                                            reduce2givenCells=T,
                                            rename_last_2cellidLettersTo = NULL, # sometimes, cell Ids in seuratobjects have different last number than cell Ids in unfiltered hd5 file. Ust this string to replace the last 2 letters in the cell IDs of the hd5 file to match cell IDs of the seuratobject
                                            beVerbose=F,
                                            testmode = F) {
  # seuratobject = mouse_filt
  # raw_h5_10x_fn = here('../../01_daten/1911_sympath_sc-seq/single_cell_berlin/results/X_WTB6J_PBS/X_WTB6J_PBS/outs/raw_feature_bc_matrix.h5')
  # experiment_name <-  "X_wt_pbs"
  # empty_drops_lower = 100
  # ribogenes = ribogenes_mouse
  # mitogenes=mitogenes_mouse


  #
  #   myrun10x = todofile_hamster$run10x[2]
  #   myzeile = todofile_hamster[run10x==myrun10x]
  #   seuratobject = hamster_filt_decontx_list[[myrun10x]]
  #   raw_h5_10x_fn = myzeile$raw_h5_10x_fn
  #   experiment_name <-  paste0("", myrun10x)
  # empty_drops_lower = 100
  # ribogenes=ribogenes_hamster
  # mitogenes=mitogenes_hamster
  # reduce2givenCells=T
  #
  # beVerbose=T
  # rename_last_2cellidLettersTo = "-2"

  # beVerbose=T
  #   myrun10x = todofile_humanCharite$run10x[1]
  #   myzeile = todofile_humanCharite[run10x==myrun10x]
  # seuratobject = humanCharite_filt_decontx_list[[myrun10x]]
  #   raw_h5_10x_fn = myzeile$raw_h5_10x_fn
  #   experiment_name <-  paste0("", myrun10x)
  # rename_last_2cellidLettersTo = "-1"
  # ribogenes=ribogenes_humanCharite
  # mitogenes=mitogenes_humanCharite
  # orthologue_genenames = orthologues_humanCharite
  # testmode =T
  library(DropletUtils)
  library(celda)

  seuratobject_name = deparse(substitute(seuratobject))

  # use DropletUtils to load data
  message(Sys.time(), '...Loading:\n', raw_h5_10x_fn )

  sce <- read10xCounts(raw_h5_10x_fn, sample.names = experiment_name, type = 'HDF5', version = 'auto')
  # label rows and columns
  colnames(sce) <- paste0(sce$Sample,"_", sce$Barcode)
  rownames(sce) <- rowData(sce)$Symbol

  if(is.null(rename_last_2cellidLettersTo)==F) colnames(sce) <- str_replace(colnames(sce), "-[0-9]$", rename_last_2cellidLettersTo)

  message(Sys.time(), "...Using experiment name '", experiment_name, "', read cell IDs of ambient RNA corrected cells will be like\n", head(colnames(sce)) %>% paste(., collapse = "\n"))
  message(Sys.time(), "...cell IDs of given Seurat object are like\n", head(colnames(seuratobject)) %>% paste(., collapse = "\n"))

  qlist64 = venn2(colnames(seuratobject), colnames(sce), plotte = testmode)
  str(qlist64)
  if(length(qlist64$q1)==0) stop("Stopping.... Please use an  'experiment_name' argument that matches cell IDs of given and read data so that cell IDs match!")

  if(testmode ==T)  sce = sce[, unique(c(colnames(sce)[1:11111], qlist64$q1))]

  # keep only droplets that have >5 expression
  message(Sys.time(), ' ...removing droplets that have less or equal 5 expression counts')
  keep_drop <- colSums(counts(sce) > 0) > 5
  sce <- sce[, keep_drop]


  # keep only genes that have nonzero expression
  message(Sys.time(), ' ...removing genes that have zero expression')

  keep_feature <- rowSums(counts(sce) > 0) > 0
  sce <- sce[keep_feature, ]
  dimension2 = dim(sce)

  message(Sys.time(), "...Filtered to  ", dimension2[1] %>% huebsch(), " genes and ", dimension2[2]%>% huebsch(), " cells")

  # split the data into cells and empties by running emptyDrops
  message(Sys.time(), '...Running emptyDrops to find cells and empties')

  # raw_sce = as.SingleCellExperiment(raw_seurat)
  empty_drops_out <- emptyDrops(assay(sce, 'counts'), lower = empty_drops_lower, ignore = 10) # If retain is not specified, it is set to the total count at the knee point detected by barcodeRanks  retain = empty_drops_retain,
  # # Users can interpret ignore as the minimum total count required for a barcode to be considered as a potential cell
  message(Sys.time(), ' ...complete')


  # get the cells and empties
  message(Sys.time(), '...Splitting dataset into cells and empties')
  is.cell <- empty_drops_out$FDR <= 0.01
  message(Sys.time(), ' ...setting NA values to FALSE')
  is.cell[is.na(is.cell)] <- FALSE
  message(Sys.time(), paste0(c(' ...found ', sum(is.cell), ' cells and ', sum(!is.cell), ' empties')))
  sce_cells <- sce[, which(is.cell)]
  sce_empty <- sce[, which(!is.cell)]
  message(Sys.time(), ' ...done')

  # run decontX with empty droplets
  message(Sys.time(), '...Running decontX')
  result <- celda::decontX(x = sce_cells, background = sce_empty, verbose = beVerbose)
  message(Sys.time(), ' ...complete')



  umap <- reducedDim(result, "decontX_UMAP")

  plot(plotDimReduceCluster(x = result$decontX_clusters,
                            dim1 = umap[, 1], dim2 = umap[, 2]) +ggtitle(seuratobject_name, subtitle = unique(seuratobject$run10x)))

  plot( plotDecontXContamination(result)  +ggtitle(seuratobject_name, subtitle = unique(seuratobject$run10x)))



  result_matrix = as.matrix(result@assays@data$decontXcounts) %>% round()

  message(Sys.time(), "...Recreating Seurat object")
  seuratobject_decontx_withbg =  CreateSeuratObject(counts = result_matrix)

  qlist65 = venn2(names(seuratobject@meta.data), names(seuratobject_decontx_withbg@meta.data), plotte = beVerbose)


  for(attribnames in setdiff(qlist65$q2,names(seuratobject_decontx_withbg@meta.data)) %>% unique()) {
    message(Sys.time(), "...trying to add original column (via barcode ID) '", attribnames, "'")
    seuratobject_decontx_withbg$addedatrib = seuratobject@meta.data[match_hk(colnames(seuratobject_decontx_withbg) , colnames(seuratobject)),attribnames]

    setnames(seuratobject_decontx_withbg@meta.data, "addedatrib", attribnames)
  }


  # add decont results
  result$decontX_umap1 = umap[, 1]
  result$decontX_umap2 = umap[, 2]
  result$id = paste0(result$Sample, "_", result$Barcode)

  toadd = c('decontX_umap1', 'decontX_umap2', 'decontX_clusters', 'decontX_contamination')

  qlist5422 = venn2(colnames(seuratobject_decontx_withbg) , result$id, plotte = beVerbose)
  # str(qlist5422)
  for(attribnames2 in toadd) {
    message(Sys.time(), "...trying to add DecontX results (via barcode ID) '", attribnames2, "'")

    seuratobject_decontx_withbg$addedatrib = result@colData[match_hk(colnames(seuratobject_decontx_withbg) , result$id),attribnames2]

    setnames(seuratobject_decontx_withbg@meta.data, "addedatrib", attribnames2)
  }

  # ggplot(seuratobject_decontx_withbg@meta.data, aes(decontX_umap1, decontX_umap2, col = decontX_contamination)) + geom_point()
  #
  # ggplot(seuratobject_decontx_withbg@meta.data, aes(decontX_umap1, decontX_umap2, col = celltype)) + geom_point()

  if(is.null(mitogenes)==F){
    message(Sys.time(), "...recalculating pct.mito")
    mitogenes_present = mitogenes[mitogenes %in% rownames(seuratobject_decontx_withbg)]
    message(Sys.time(), "...using for recalculation pct.mito genes ", paste(sort(mitogenes_present), collapse = ", "))

    if(is.null(orthologue_genenames) ==T) seuratobject_decontx_withbg$pct.mito <- PercentageFeatureSet(seuratobject_decontx_withbg,features = mitogenes_present) else seuratobject_decontx_withbg$pct.mito <- PercentageFeatureSet(seuratobject_decontx_withbg[rownames(seuratobject_decontx_withbg) %in% orthologue_genenames,],features = mitogenes_present)
  }

  if(is.null(ribogenes)==F) {
    message(Sys.time(), "...recalculating pct.ribo")
    ribogenes_present = ribogenes[ribogenes %in% rownames(seuratobject_decontx_withbg)]
    message(Sys.time(), "...using for recalculation pct.mito genes ", paste(sort(ribogenes_present), collapse = ", "))

    if(is.null(orthologue_genenames) ==T) seuratobject_decontx_withbg$pct.ribo =  PercentageFeatureSet(seuratobject_decontx_withbg, features = ribogenes_present) else seuratobject_decontx_withbg$pct.ribo =  PercentageFeatureSet(seuratobject_decontx_withbg[rownames(seuratobject_decontx_withbg) %in% orthologue_genenames,], features = ribogenes_present)
  }

  message(Sys.time(), "...trying to match original column (via barcode ID) '", attribnames, "'")


  qlist66 = venn2(colnames(seuratobject), colnames(seuratobject_decontx_withbg) , mylabels = c(seuratobject_name, "Ambient-RNA-removed object"), mytitle = paste0("\nNumber of Dropletts with Cells\n", unique(result$Sample)))

  if(reduce2givenCells==T) {
    message(Sys.time(), "...Reducing to cells of the initial object....")
    seuratobject_decontx_withbg = seuratobject_decontx_withbg[, qlist66$q1]
  }

  seuratobject_decontx_withbg
}

doSCTransform = function(seuratobject, experimentIDcol = "run10x", variable.features.n = 4000, ...) {
  # seuratobject = mouse;experimentIDcol = "run10x"

  message(Sys.time(), "...Running Seurat SCT separately on following cells:")
  mytable(seuratobject[[experimentIDcol]])

  n_experiments = uniqueN(names(table(seuratobject[[experimentIDcol]])))

  if(n_experiments ==1) seuratobject_list = list(seuratobject) else seuratobject_list <- SplitObject(object = seuratobject, split.by = experimentIDcol)
  seuratobject_list


  for (i in 1:length(seuratobject_list)) {
    # i=1
    message(Sys.time(), "...SCT on ", i, "\n-------------------------------------------------------")
    seuratobject_list[[i]] <- SCTransform(seuratobject_list[[i]], verbose = T,variable.features.n = variable.features.n, ... ) # TODO include "G2M.Score", "S.Score" regression if necessary
  }
  print(table(names(warnings() ))) #
  message(Sys.time(), '...Notes about warning:\n"iteration limit reached"\nChristophH commented on 22 May 2019 - These warnings are showing that there are some genes for which it is hard to reliably estimate theta (presumably because of very few non-zero observations). Usually we don´t worry about these warnings too much, since we regularize the parameters in a later step, thus averaging out uncertainty of individual gene parameters. https://github.com/ChristophH/sctransform/issues/25')


  list.features <- SelectIntegrationFeatures(object.list = seuratobject_list, nfeatures = 4000)

  if(length(seuratobject_list)>1) {
    seuratobject <- merge(seuratobject_list[[1]],
                          y = seuratobject_list[2:length(seuratobject_list)],
                          # project = "seuratobject",
                          merge.data = TRUE)

  } else seuratobject = seuratobject_list[[1]]

  VariableFeatures(seuratobject) <- list.features
  seuratobject
}

plot3clusterings = function(seuratobject, clustervar1 = 'SCT_snn_res.0.2', clustervar2 = 'SCT_snn_res.0.4', clustervar3 = 'SCT_snn_res.0.8') {
  # clustervar1 = 'SCT_snn_res.0.2'; clustervar2 = 'SCT_snn_res.0.4'; clustervar3 = 'SCT_snn_res.0.8'
  seuratobject_name = deparse(substitute(seuratobject))
  p1 = DimPlot(seuratobject, group.by = clustervar1, label = T)
  p2 = DimPlot(seuratobject, group.by = clustervar2, label = T)
  p3 = DimPlot(seuratobject, group.by = clustervar3, label = T)

  plot(p1+p2+p3 + patchwork::plot_annotation(title = seuratobject_name))
}


calcDoubletts <- function(seuratobject,
                          doSCT = T,  # was the seurat object preprocessed with SCT, see # maximal number of Pcs, see https://github.com/chris-mcginnis-ucsf/DoubletFinder
                          maxPC = 10, # maximal number of Pcs, see https://github.com/chris-mcginnis-ucsf/DoubletFinder
                          annotation_column = "SCT_snn_res.0.4" , # Clusternames

                          doublets_expected_fromLoading = NULL, # # corresponding expectation if 8000 cells loaded, see  https://assets.ctfassets.net/an68im79xiti/4tjk4KvXzTWgTs8f3tvUjq/2259891d68c53693e753e1b45e42de2d/CG000183_ChromiumSingleCell3__v3_UG_Rev_C.pdf
                          seed = 2712, # seed for reproduction
                          experimental_batchID = "run10x"
) {
  # require(jtools)
  seuratobject_name = deparse(substitute(seuratobject_name))
  library(DoubletFinder) # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

  # seuratobject = monkey_filt_decontx
  # annotation_column = "SCT_snn_res.0.2" # reasonable number of cluster, should roughly match number of cell types
  # doublets_expected_fromLoading = 0.03
  #   experimental_batchID = "run10x"
  set.seed(seed)
  stopifnot(experimental_batchID %in% names(seuratobject@meta.data))
  message(Sys.time(), "...expecting '",experimental_batchID,"' as sub-experiment-identifyer...")
  mytable(seuratobject@meta.data[[experimental_batchID]])
  stopifnot(annotation_column %in% names(seuratobject@meta.data))
  message(Sys.time(), "...using column '",annotation_column,"' for clustering...")
  orig.idents = unique(seuratobject@meta.data[[experimental_batchID]])
  message(paste0(Sys.time(), "...Running DoubletFinder for\n", paste(orig.idents, collapse = "\n")))
  # separate for each

  if("pANN" %in% names(seuratobject@meta.data)) seuratobject$pANN = NULL

  totalres = lapply(orig.idents, function(myident) {
    # myident = orig.idents[1]
    res = c()
    message(Sys.time(), "...=======================================================\nWorking on ", myident)

    subseurat  = seuratobject[,seuratobject@meta.data[[experimental_batchID]]== myident]
    subseurat


    # https://github.com/chris-mcginnis-ucsf/DoubletFinder

    ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
    sweep.res.list_subseurat <- paramSweep_v3(subseurat, PCs = 1:maxPC, sct = doSCT)
    sweep.stats_subseurat <- summarizeSweep(sweep.res.list_subseurat, GT = FALSE)
    bcmvn_subseurat <- find.pK(sweep.stats_subseurat)
    bcmvn_subseurat
    # plot(bcmvn_subseurat)
    mypk = bcmvn_subseurat[bcmvn_subseurat$BCmetric==max(bcmvn_subseurat$BCmetric), "pK"] %>% as.character %>% as.numeric
    mypk
    message("Using pk of ",mypk)
    res$bcmvn = bcmvn_subseurat
    res$pk =mypk

    ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
    # sweep.res.list_subseurat <- paramSweep_v3(subseurat, PCs = 1:maxPC, sct = FALSE)
    # gt.calls <- subseurat@meta.data[rownames(sweep.res.list_subseurat[[1]]), "GT"].   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results
    # sweep.stats_subseurat <- summarizeSweep(sweep.res.list_subseurat, GT = TRUE, GT.calls = gt.calls)
    # bcmvn_subseurat <- find.pK(sweep.stats_subseurat)

    ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
    annotations <- subseurat@meta.data[,annotation_column]
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- subseurat@meta.data$ClusteringResults
    homotypic.prop


    if(is.null(doublets_expected_fromLoading)) {
      cellspresent =  dim(subseurat)[2]
      estimated_doublettrate = 0.0007665152*cellspresent
      message(Sys.time(), "...estimating Doublett rate from number of present cells per batch ", "(",cellspresent, " -> ca. ", round(estimated_doublettrate, 1), "%)")
      doublettrate = estimated_doublettrate /100
    } else doublettrate = doublets_expected_fromLoading
    message(Sys.time(), "...Using total doublet rate of ",doublettrate, " (For 10x data, should reflect loading acc. to https://assets.ctfassets.net/an68im79xiti/4tjk4KvXzTWgTs8f3tvUjq/2259891d68c53693e753e1b45e42de2d/CG000183_ChromiumSingleCell3__v3_UG_Rev_C.pdf" )

    nExp_poi <- round(doublettrate*nrow(subseurat@meta.data))

    nExp_poi
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    nExp_poi.adj

    res$nExp_poi = nExp_poi
    res$nExp_poi.adj = nExp_poi.adj

    message(Sys.time(), "...Running DoubletFinder with varying classification stringencies ----------------------------------------------------------------")
    subseurat <- doubletFinder_v3(subseurat, PCs = 1:maxPC, pN = 0.25, pK = mypk, nExp = nExp_poi, reuse.pANN = FALSE, sct = doSCT)
    pann_name = grep("^pANN", names(subseurat@meta.data), value = T)
    pann_name
    res$pann_name = pann_name
    # p_pANN = Seurat::FeaturePlot(subseurat, features =  pann_name) +
    #   ggtitle(myident, subtitle = pann_name)
    # plot(p_pANN)
    # res$p_pANN = p_pANN

    subseurat <- doubletFinder_v3(subseurat, PCs = 1:maxPC, pN = 0.25, pK = mypk, nExp = nExp_poi.adj, reuse.pANN = pann_name, sct = doSCT)

    subseurat_attrib = subseurat@meta.data %>% as.data.table(.,keep.rownames = T)

    DF_finalnames = grep("DF\\.classific", names(subseurat@meta.data), value = T)
    DF_finalnames


    # p_notConsidHomoDoubl = Seurat::DimPlot(subseurat, group.by =   DF_finalnames[1]) +
    #   ggtitle(myident, subtitle = DF_finalnames[1])
    # p_ConsidHomoDoubl = Seurat::DimPlot(subseurat, group.by =   DF_finalnames[2])+
    #   ggtitle(myident, subtitle = DF_finalnames[2])
    #
    # res$p_notConsidHomoDoubl = p_notConsidHomoDoubl
    # res$p_ConsidHomoDoubl = p_ConsidHomoDoubl


    subseurat_attrib[, DF.classifications_notConsidHomoDoubl := get(DF_finalnames[1])]
    subseurat_attrib[, DF.classifications_ConsidHomoDoubl := get(DF_finalnames[2])]

    subseurat_attrib[,.N, .(DF.classifications_notConsidHomoDoubl, DF.classifications_ConsidHomoDoubl)]
    subseurat_attrib[DF.classifications_notConsidHomoDoubl != DF.classifications_ConsidHomoDoubl,.N , get(annotation_column)]

    annodata_pre =  subseurat_attrib[,.N, .(get(annotation_column),DF.classifications_ConsidHomoDoubl)]
    names(annodata_pre)[names(annodata_pre)=="get"] = annotation_column
    annodata  = dcast.data.table(annodata_pre, get(annotation_column) ~ DF.classifications_ConsidHomoDoubl, value.var = "N", fill = 0)
    names(annodata)[names(annodata)=="annotation_column"] = annotation_column
    annodata[,proz_doublet := paste0((Doublet/(Doublet + Singlet)) %>% proz(.,stellen = 0), "\n(",Doublet, ")")]
    annodata2 = cbind(DF.classifications_ConsidHomoDoubl = "Doublet", annodata)

    p_bar = ggplot(subseurat_attrib, aes_string(x = annotation_column, fill = 'DF.classifications_ConsidHomoDoubl')) +
      geom_bar(position = "fill", alpha = 0.8) + scale_y_continuous(labels = label_percent(accuracy = 1), breaks = (0:10)/10) + theme_hc() +
      theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1),
            legend.position = "top",
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
      ) +
      scale_fill_colorblind()+
      labs(fill = "")+
      geom_text(data =annodata2, aes_string(x = annotation_column, label = 'proz_doublet', y = 1.1), col = viridis_pal()(1), size =3) +
      ggtitle(paste(seuratobject_name , "-", myident), subtitle = 'DF.classifications_ConsidHomoDoubl')

    plot(p_bar)

    res$p_barplot = p_bar
    res$annotation= subseurat_attrib
    message(Sys.time(), "...DONE Doublet analysis for ", myident)
    res


  }
  )
  message(Sys.time(), "...Collecting Doublet results back into seurat object... ")
  names(totalres) = orig.idents


  doubletdf = lapply(totalres, function(x) {
    setnames(x$annotation, grep("pANN", names(x$annotation), value = T), "pANN")
  }) %>% rbindlist(., fill = T) %>% .[,grep('DF\\.classifications_0', names(.), invert = T), with = F]

  doubletdf

  seuratobject$pANN = doubletdf[toolboxH::match_hk(colnames(seuratobject), doubletdf$rn),pANN]
  seuratobject$DF.classifications_notConsidHomoDoubl = doubletdf[toolboxH::match_hk(colnames(seuratobject), doubletdf$rn),DF.classifications_notConsidHomoDoubl]
  seuratobject$DF.classifications_ConsidHomoDoubl = doubletdf[toolboxH::match_hk(colnames(seuratobject), doubletdf$rn),DF.classifications_ConsidHomoDoubl]



  plot(DimPlot(seuratobject, group.by = 'DF.classifications_ConsidHomoDoubl'))

  seuratobject

}

correctAmbientRNA_noEmptyCells = function(seuratobject,  ribogenes=NULL, mitogenes=NULL ,
                                          orthologue_genenames = NULL, # if mitogenes or ribogenes are provided than  if also providing orthologue_genenames, percentage of ribogenes and mitogenes relates to orthologue_genenames also included in the seuratobject, only
                                          beVerbose = F ) {
  # seuratobject = monkey_filt

  # ribogenes = ribogenes_monkey
  # beVerbose = T
  # mitogenes=mitogenes_monkey
  #   orthologue_genenames = orthologues_monkey
  # seuratobject = SplitObject(object = monkey_filt, split.by = "run10x")[[1]]


  seuratobject_name = deparse(substitute(seuratobject))

  require(celda)
  require(SingleCellExperiment)
  seuratobject_sce = as.SingleCellExperiment(seuratobject)

  # run decontX with empty droplets
  message(Sys.time(), '...Running decontX')

  seuratobject_sce <- decontX(seuratobject_sce, verbose = beVerbose)


  message(Sys.time(), ' ...complete')
  umap <- reducedDim(seuratobject_sce, "decontX_UMAP")

  plot(plotDimReduceCluster(x = seuratobject_sce$decontX_clusters,
                            dim1 = umap[, 1], dim2 = umap[, 2]) +ggtitle(seuratobject_name, subtitle = unique(seuratobject$run10x)))

  plot( plotDecontXContamination(seuratobject_sce)  +ggtitle(seuratobject_name, subtitle = unique(seuratobject$run10x)))

  message(Sys.time(), "...Recreating Seurat object")
  seuratobject_decontx_nobg =  CreateSeuratObject(counts = seuratobject_sce@assays@data$decontXcounts %>% round())

  qlist65 = venn2(names(seuratobject@meta.data), names(seuratobject_decontx_nobg@meta.data), plotte = F)

  for(attribnames in qlist65$q2) {
    # attribnames = qlist65$q2[1]
    message(Sys.time(), "...trying to add original column (via barcode ID) '", attribnames, "'")
    seuratobject_decontx_nobg$addedatrib = seuratobject@meta.data[match_hk(colnames(seuratobject_decontx_nobg) , colnames(seuratobject)),attribnames]
    setnames(seuratobject_decontx_nobg@meta.data, "addedatrib", attribnames)
  }


  ##################################
  # add decont results
  seuratobject_sce$decontX_umap1 = umap[, 1]
  seuratobject_sce$decontX_umap2 = umap[, 2]

  toadd = c('decontX_umap1', 'decontX_umap2', 'decontX_clusters', 'decontX_contamination')


  for(attribnames2 in toadd) {
    message(Sys.time(), "...trying to add DecontX results (via barcode ID) '", attribnames2, "'")

    seuratobject_decontx_nobg$addedatrib = seuratobject_sce@colData[match_hk(colnames(seuratobject_decontx_nobg) , colnames(seuratobject_sce)),attribnames2]

    setnames(seuratobject_decontx_nobg@meta.data, "addedatrib", attribnames2)
  }



  if(is.null(mitogenes)==F){
    message(Sys.time(), "...recalculating pct.mito")
    mitogenes_present = mitogenes[mitogenes %in% rownames(seuratobject_decontx_nobg)]
    message(Sys.time(), "...using for recalculation pct.mito genes ", paste(sort(mitogenes_present), collapse = ", "))

    if(is.null(orthologue_genenames) ==T) seuratobject_decontx_nobg$pct.mito <- PercentageFeatureSet(seuratobject_decontx_nobg,features = mitogenes_present) else seuratobject_decontx_nobg$pct.mito <- PercentageFeatureSet(seuratobject_decontx_nobg[rownames(seuratobject_decontx_nobg) %in% orthologue_genenames,],features = mitogenes_present)
  }

  if(is.null(ribogenes)==F) {
    message(Sys.time(), "...recalculating pct.ribo")
    ribogenes_present = ribogenes[ribogenes %in% rownames(seuratobject_decontx_nobg)]
    message(Sys.time(), "...using for recalculation pct.mito genes ", paste(sort(ribogenes_present), collapse = ", "))

    if(is.null(orthologue_genenames) ==T) seuratobject_decontx_nobg$pct.ribo =  PercentageFeatureSet(seuratobject_decontx_nobg, features = ribogenes_present) else seuratobject_decontx_nobg$pct.ribo =  PercentageFeatureSet(seuratobject_decontx_nobg[rownames(seuratobject_decontx_nobg) %in% orthologue_genenames,], features = ribogenes_present)
  }


  seuratobject_decontx_nobg
}



# # >>MOUSE----

mouse = readRDS(here("data/s511_1_diet_mouse.RDS"))
mouse

mouse = DietSeurat(mouse)
# ensuring basic QC ----
qlist_mouse = venn3(rownames(mouse), orthologues$`Mouse gene name`, orthologues$`Mouse gene stable ID`)


mitogenes_mouse = orthologues[grep('^MT-', `Human gene name`), unique(c(`Mouse gene name`,`Mouse gene stable ID`))] %>% na.omit() %>% intersect(rownames(mouse), .)
mitogenes_mouse

ribogenes_mouse = orthologues[grep('^RP[SL][0-9]*$', `Human gene name`), unique(c(`Mouse gene name`,`Mouse gene stable ID`))] %>% na.omit() %>% intersect(rownames(mouse), .)
ribogenes_mouse

s.genes_mouse <- orthologues[`Human gene name` %in% cc.genes$s.genes ,unique(c(`Mouse gene name`,`Mouse gene stable ID`))] %>% na.omit() %>% intersect(rownames(mouse), .)
s.genes_mouse

g2m.genes_mouse <- orthologues[`Human gene name` %in% cc.genes$g2m.genes,unique(c(`Mouse gene name`,`Mouse gene stable ID`))] %>% na.omit() %>% intersect(rownames(mouse), .)
g2m.genes_mouse


orthologues_mouse = unique(na.omit(c(orthologues$`Mouse gene name`, orthologues$`Mouse gene stable ID`)))

mouse$pct.mito <- PercentageFeatureSet(mouse[rownames(mouse) %in% orthologues_mouse,],features = mitogenes_mouse) #

mouse$pct.ribo =  PercentageFeatureSet(mouse[rownames(mouse) %in% orthologues_mouse,], features = ribogenes_mouse)

mouse@meta.data$run10x = paste0("mouse_",mouse@meta.data$orig.ident)
table(mouse@meta.data$run10x)
mouse@meta.data$species = "mouse"

# here's some basic statistics.
stats1_mouse  = mouse@meta.data %>%
  dplyr::group_by(run10x) %>%
  dplyr::summarise(ncells=dplyr::n(),
                   # ndoublets=sum(scrublet_prediction=='doublet'),
                   nGene=median(nFeature_RNA),
                   nUMI=median(nCount_RNA),
                   pct.mito=mean(pct.mito),
                   pct.ribo=mean(pct.ribo))

stats1_mouse


plotQC(mouse)




# Basic QC filter


mouse_filt = doBAsicFiltering(mouse)

# Ambient RNA background correction with available emptydroplet data----






mouse_filt_decontx = correctAmbientRNA_inclEmptyCells(seuratobject=mouse_filt,
                                                      raw_h5_10x_fn = here('data/mouse_seurat_raw/raw_feature_bc_matrix.h5'),
                                                      experiment_name <-  "X_wt_pbs",
                                                      empty_drops_lower = 100,
                                                      ribogenes=ribogenes_mouse,
                                                      mitogenes=mitogenes_mouse ,
                                                      orthologue_genenames = orthologues_mouse,
                                                      reduce2givenCells=T)

# Clustern, to identifyt  Doubletts ----
# refilter as ambient RNA correction might have changed data



mouse_filt_decontx2 <- subset(mouse_filt_decontx, (pct.mito < param_QC_filter_mito_max) &
                                (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                                (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                                (nCount_RNA < param_QC_filter_RNAcounts_max) &
                                (nCount_RNA > param_QC_filter_RNAcounts_min) ) %>% # & (scrublet_prediction!='doublet')
  doSCTransform(seuratobject = .,experimentIDcol = "run10x", variable.features.n = 4000,   vars.to.regress = "pct.mito") %>%
  CellCycleScoring(s.features = s.genes_mouse, g2m.features = g2m.genes_mouse,  assay = "SCT") %>%
  RunPCA(features=VariableFeatures(.), assay = "SCT", npcs = 50,verbose=TRUE)%>%
  RunUMAP(dims=1:30,verbose=TRUE,  assay = "SCT" ) %>%
  FindNeighbors(dims=1:30, assay = "SCT", verbose=TRUE)  %>%
  FindClusters( resolution = 0.8) %>%
  FindClusters( resolution = 0.2) %>%
  FindClusters( resolution = 0.4)
mouse_filt_decontx2

# plot Cluster ---


plot3clusterings(mouse_filt_decontx2)


mouse_filt_decontx2 = calcDoubletts(seuratobject = mouse_filt_decontx2,
                                    annotation_column = "SCT_snn_res.0.4", # reasonable number of cluster, should roughly match number of cell types
                                    doublets_expected_fromLoading = 0.039 # corresponding expectation if 8000 cells loaded, see  https://assets.ctfassets.net/an68im79xiti/4tjk4KvXzTWgTs8f3tvUjq/2259891d68c53693e753e1b45e42de2d/CG000183_ChromiumSingleCell3__v3_UG_Rev_C.pdf
)
mouse_filt_decontx2

# Filter again including doublets


mouse_filt_decontx3 <- subset(mouse_filt_decontx2, (pct.mito < param_QC_filter_mito_max) &
                                (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                                (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                                (nCount_RNA < param_QC_filter_RNAcounts_max) &
                                (nCount_RNA > param_QC_filter_RNAcounts_min) &
                                DF.classifications_ConsidHomoDoubl == "Singlet")


mouse_filt_decontx3

saveRDS(mouse_filt_decontx2, file = here('results/s602_9v2_mouse_preprocessed_withDoubletts.RDS'))

# >>MONKEY----

monkey = readRDS(here("data/s101_2_greenmonkey_GSM4743527_16434_16436.RDS"))
monkey

monkey = DietSeurat(monkey)

# replace "" with NA in Gene names----
showNA(orthologues, showAllNoNA = F)
for(i in c( "Vervet-AGM gene name", 'Vervet-AGM gene stable ID')) {
  orthologues[get(i)=="",(i):=NA]
}
showNA(orthologues, showAllNoNA = F)


# ensuring basic QC ----
qlist_monkey1 = venn3(rownames(monkey), orthologues$`Vervet-AGM gene name`, orthologues$`Vervet-AGM gene stable ID`)

mitogenes_monkey = orthologues[grep('^MT-', `Human gene name`), unique(c(`Vervet-AGM gene stable ID`, `Vervet-AGM gene name`))] %>% na.omit() %>% intersect(rownames(monkey), .)
mitogenes_monkey

ribogenes_monkey = orthologues[grep('^RP[SL][0-9]*$', `Human gene name`), unique(c(`Vervet-AGM gene stable ID`, `Vervet-AGM gene name`))] %>% na.omit() %>% intersect(rownames(monkey), .)
ribogenes_monkey

s.genes_monkey <- orthologues[`Human gene name` %in% cc.genes$s.genes ,  unique(c(`Vervet-AGM gene stable ID`, `Vervet-AGM gene name`))]  %>% na.omit() %>% intersect(rownames(monkey), .)
s.genes_monkey

g2m.genes_monkey <- orthologues[`Human gene name` %in% cc.genes$g2m.genes, unique(c(`Vervet-AGM gene stable ID`, `Vervet-AGM gene name`))]  %>% na.omit() %>% intersect(rownames(monkey), .)
g2m.genes_monkey


orthologues_monkey = unique(na.omit(c(orthologues$`Vervet-AGM gene name`, orthologues$`Vervet-AGM gene stable ID`)))


monkey$pct.mito <- PercentageFeatureSet(monkey[rownames(monkey) %in% orthologues_monkey,],features = mitogenes_monkey)
monkey$pct.ribo =  PercentageFeatureSet(monkey[rownames(monkey) %in% orthologues_monkey,], features = ribogenes_monkey)

monkey@meta.data$run10x = paste0("monkey_",monkey@meta.data$orig.ident)
table(monkey@meta.data$run10x)
monkey@meta.data$species = "monkey"

# here's some basic statistics.
stats1_monkey  = monkey@meta.data %>%
  dplyr::group_by(run10x) %>%
  dplyr::summarise(ncells=dplyr::n(),
                   # ndoublets=sum(scrublet_prediction=='doublet'),
                   nGene=median(nFeature_RNA),
                   nUMI=median(nCount_RNA),
                   pct.mito=mean(pct.mito),
                   pct.ribo=mean(pct.ribo))

stats1_monkey
plotQC(monkey)



# Basic QC filter
monkey_filt = doBAsicFiltering(monkey)

# correct for contamination with ambient RNA ----

monkey_filt_decontx_list = SplitObject(object = monkey_filt, split.by = "run10x") %>%   lapply(., function(x) correctAmbientRNA_noEmptyCells(seuratobject = x,
                                                                                                                                             ribogenes = ribogenes_monkey,
                                                                                                                                             mitogenes = mitogenes_monkey,
                                                                                                                                             orthologue_genenames = orthologues_monkey))

monkey_filt_decontx = merge(monkey_filt_decontx_list[[1]],y = monkey_filt_decontx_list[2:length(monkey_filt_decontx_list)], merge.data = TRUE)
monkey_filt_decontx

# Identify  Doubletts ----
# refilter as ambient RNA correction might have changed data


monkey_filt_decontx2 <- subset(monkey_filt_decontx, (pct.mito < param_QC_filter_mito_max) &
                                 (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                                 (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                                 (nCount_RNA < param_QC_filter_RNAcounts_max) &
                                 (nCount_RNA > param_QC_filter_RNAcounts_min) ) %>% # & (scrublet_prediction!='doublet')
  doSCTransform(seuratobject = .,experimentIDcol = "run10x", variable.features.n = 4000,   vars.to.regress = "pct.mito")  %>%
  CellCycleScoring(s.features = s.genes_monkey, g2m.features = g2m.genes_monkey,  assay = "SCT") %>%
  RunPCA(features=VariableFeatures(.), assay = "SCT", npcs = 50,verbose=TRUE)%>%
  RunUMAP(dims=1:30,verbose=TRUE,  assay = "SCT" ) %>%
  FindNeighbors(dims=1:30, assay = "SCT", verbose=TRUE)  %>%
  FindClusters( resolution = 0.8) %>%
  FindClusters( resolution = 0.2) %>%
  FindClusters( resolution = 0.4)

monkey_filt_decontx2
# plot Cluster ---
plot3clusterings(monkey_filt_decontx2)



monkey_filt_decontx2 = calcDoubletts(seuratobject = monkey_filt_decontx2,
                                     annotation_column = "SCT_snn_res.0.2", # reasonable number of cluster, should roughly match number of cell types
                                     doublets_expected_fromLoading = NULL # corresponding expectation  see  https://assets.ctfassets.net/an68im79xiti/4tjk4KvXzTWgTs8f3tvUjq/2259891d68c53693e753e1b45e42de2d/CG000183_ChromiumSingleCell3__v3_UG_Rev_C.pdf page 17
)
monkey_filt_decontx2

# Filter again including doublets

monkey_filt_decontx3 <- subset(monkey_filt_decontx2, (pct.mito < param_QC_filter_mito_max) &
                                 (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                                 (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                                 (nCount_RNA < param_QC_filter_RNAcounts_max) &
                                 (nCount_RNA > param_QC_filter_RNAcounts_min) &
                                 DF.classifications_ConsidHomoDoubl == "Singlet")

monkey_filt_decontx3


saveRDS(monkey_filt_decontx2, file = here('results/s602_9v2_monkey_preprocessed_withDoubletts.RDS'))

# # >>PIG----
pig = readRDS(here("data/s101_2_pig_DS026_DS027.RDS"))
pig

pig = DietSeurat(pig)

# replace "" with NA in Gene names----
showNA(orthologues, showAllNoNA = F)
for(i in c( "Pig gene name", 'Pig gene stable ID')) {
  orthologues[get(i)=="",(i):=NA]
}
showNA(orthologues, showAllNoNA = F)


# ensuring basic QC ----
qlist_pig1 = venn3(rownames(pig), orthologues$`Pig gene name`, orthologues$`Pig gene stable ID`)

mitogenes_pig = orthologues[grep('^MT-', `Human gene name`), unique(c(`Pig gene stable ID`, `Pig gene name`))] %>% na.omit() %>% intersect(rownames(pig), .)
mitogenes_pig

ribogenes_pig = orthologues[grep('^RP[SL][0-9]*$', `Human gene name`), unique(c(`Pig gene stable ID`, `Pig gene name`))] %>% na.omit() %>% intersect(rownames(pig), .)
ribogenes_pig

s.genes_pig <- orthologues[`Human gene name` %in% cc.genes$s.genes ,  unique(c(`Pig gene stable ID`, `Pig gene name`))]  %>% na.omit() %>% intersect(rownames(pig), .)
s.genes_pig

g2m.genes_pig <- orthologues[`Human gene name` %in% cc.genes$g2m.genes, unique(c(`Pig gene stable ID`, `Pig gene name`))]  %>% na.omit() %>% intersect(rownames(pig), .)
g2m.genes_pig


orthologues_pig = unique(na.omit(c(orthologues$`Pig gene name`, orthologues$`Pig gene stable ID`)))


pig$pct.mito <- PercentageFeatureSet(pig[rownames(pig) %in% orthologues_pig,],features = mitogenes_pig)
pig$pct.ribo =  PercentageFeatureSet(pig[rownames(pig) %in% orthologues_pig,], features = ribogenes_pig)

pig@meta.data$run10x = paste0("pig_",pig@meta.data$orig.ident)
table(pig@meta.data$run10x)
pig@meta.data$species = "pig"

# here's some basic statistics.
stats1_pig  = pig@meta.data %>%
  dplyr::group_by(run10x) %>%
  dplyr::summarise(ncells=dplyr::n(),
                   # ndoublets=sum(scrublet_prediction=='doublet'),
                   nGene=median(nFeature_RNA),
                   nUMI=median(nCount_RNA),
                   pct.mito=mean(pct.mito),
                   pct.ribo=mean(pct.ribo))

stats1_pig
plotQC(pig)



# Basic QC filter
pig_filt = doBAsicFiltering(pig)

# correct for contamination with ambient RNA ----


pig_filt_decontx_list = SplitObject(object = pig_filt, split.by = "run10x") %>%   lapply(., function(x) correctAmbientRNA_noEmptyCells(seuratobject = x,
                                                                                                                                       ribogenes = ribogenes_pig,
                                                                                                                                       mitogenes = mitogenes_pig,
))

pig_filt_decontx = merge(pig_filt_decontx_list[[1]],y = pig_filt_decontx_list[2:length(pig_filt_decontx_list)], merge.data = TRUE)
pig_filt_decontx

# Identify  Doubletts ----
# refilter as ambient RNA correction might have changed data


pig_filt_decontx2 <- subset(pig_filt_decontx, (pct.mito < param_QC_filter_mito_max) &
                              (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                              (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                              (nCount_RNA < param_QC_filter_RNAcounts_max) &
                              (nCount_RNA > param_QC_filter_RNAcounts_min) ) %>% # & (scrublet_prediction!='doublet')
  doSCTransform(seuratobject = .,experimentIDcol = "run10x", variable.features.n = 4000,   vars.to.regress = "pct.mito")  %>%
  CellCycleScoring(s.features = s.genes_pig, g2m.features = g2m.genes_pig,  assay = "SCT") %>%
  RunPCA(features=VariableFeatures(.), assay = "SCT", npcs = 50,verbose=TRUE)%>%
  RunUMAP(dims=1:30,verbose=TRUE,  assay = "SCT" ) %>%
  FindNeighbors(dims=1:30, assay = "SCT", verbose=TRUE)  %>%
  FindClusters( resolution = 0.8) %>%
  FindClusters( resolution = 0.2) %>%
  FindClusters( resolution = 0.4)

pig_filt_decontx2

# plot Cluster ---
plot3clusterings(pig_filt_decontx2)



pig_filt_decontx2 = calcDoubletts(seuratobject = pig_filt_decontx2,
                                  annotation_column = "SCT_snn_res.0.2", # reasonable number of cluster, should roughly match number of cell types
                                  doublets_expected_fromLoading = NULL # corresponding expectation if 8000 cells loaded, see  https://assets.ctfassets.net/an68im79xiti/4tjk4KvXzTWgTs8f3tvUjq/2259891d68c53693e753e1b45e42de2d/CG000183_ChromiumSingleCell3__v3_UG_Rev_C.pdf page 17
)
pig_filt_decontx2

# Filter again including doublets

pig_filt_decontx3 <- subset(pig_filt_decontx2, (pct.mito < param_QC_filter_mito_max) &
                              (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                              (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                              (nCount_RNA < param_QC_filter_RNAcounts_max) &
                              (nCount_RNA > param_QC_filter_RNAcounts_min) &
                              DF.classifications_ConsidHomoDoubl == "Singlet")
pig_filt_decontx3


saveRDS(pig_filt_decontx2, file = here('results/s602_9v2_pig_preprocessed_withDoubletts.RDS'))

# >>RAT ----
rat = readRDS(file = here("data/s101_2_rat_DS024_DS025.RDS"))
rat
rat = DietSeurat(rat)

# replace "" with NA in Gene names----
showNA(orthologues, showAllNoNA = F)
for(i in c( "Rat gene name", "Rat gene stable ID")) {
  orthologues[get(i)=="",(i):=NA]
}
showNA(orthologues, showAllNoNA = F)


# ensuring basic QC ----
qlist_rat0 = venn3(rownames(rat), orthologues$`Rat gene name`, orthologues$`Rat gene stable ID`)
str(qlist_rat0)

qlist_rat1 = venn3(rownames(rat)%>% toupper(), orthologues$`Rat gene name`%>% toupper(), orthologues$`Rat gene stable ID`%>% toupper())
str(qlist_rat1)
table(str_sub(orthologues$`Rat gene stable ID`, 1,6))

# change gene names to upper case to overcome discrepancies with orhtologue names
orthologues$`Rat gene name uc` = toupper(orthologues$`Rat gene name` )

rat = Seurat.utils::RenameGenesSeurat(obj = rat, newnames = toupper(rownames(rat))) # # from https://github.com/vertesy/Seurat.utils/
rat

qlist_rat02 = venn3(rownames(rat), orthologues$`Rat gene name uc`, orthologues$`Rat gene stable ID`)
str(qlist_rat02)

mitogenes_rat = orthologues[grep('^MT-', `Human gene name`), unique(c(`Rat gene stable ID`, `Rat gene name uc`)) %>% toupper()] %>% na.omit() %>% intersect(rownames(rat)%>% toupper(), .)
mitogenes_rat = rownames(rat)[rownames(rat)%>% toupper() %in% mitogenes_rat]
mitogenes_rat

ribogenes_rat = orthologues[grep('^RP[SL][0-9]*$', `Human gene name`), unique(c(`Rat gene stable ID`, `Rat gene name uc`))%>% toupper()] %>% na.omit() %>% intersect(rownames(rat)%>% toupper(), .)
ribogenes_rat = rownames(rat)[rownames(rat)%>% toupper() %in% ribogenes_rat]
ribogenes_rat

s.genes_rat <- orthologues[`Human gene name` %in% cc.genes$s.genes ,  unique(c(`Rat gene stable ID`, `Rat gene name uc`))%>% toupper()]  %>% na.omit() %>% intersect(rownames(rat)%>% toupper(), .)
s.genes_rat = rownames(rat)[rownames(rat)%>% toupper() %in% s.genes_rat]
s.genes_rat

g2m.genes_rat <- orthologues[`Human gene name` %in% cc.genes$g2m.genes, unique(c(`Rat gene stable ID`, `Rat gene name uc`))%>% toupper()]  %>% na.omit() %>% intersect(rownames(rat)%>% toupper(), .)
g2m.genes_rat = rownames(rat)[rownames(rat)%>% toupper() %in% g2m.genes_rat]
g2m.genes_rat


orthologues_rat = unique(na.omit(c(orthologues$`Rat gene name uc`, orthologues$`Rat gene stable ID`)))

mitogenes_rat %in% rownames(rat[rownames(rat) %in% orthologues_rat,])

rat$pct.mito <- PercentageFeatureSet(rat[rownames(rat) %in% orthologues_rat,],features = mitogenes_rat)
rat$pct.ribo =  PercentageFeatureSet(rat[rownames(rat) %in% orthologues_rat,], features = ribogenes_rat)

rat@meta.data$run10x = paste0("rat_",rat@meta.data$orig.ident)
table(rat@meta.data$run10x)
rat@meta.data$species = "rat"

# here's some basic statistics.
stats1_rat  = rat@meta.data %>%
  dplyr::group_by(run10x) %>%
  dplyr::summarise(ncells=dplyr::n(),
                   # ndoublets=sum(scrublet_prediction=='doublet'),
                   nGene=median(nFeature_RNA),
                   nUMI=median(nCount_RNA),
                   pct.mito=mean(pct.mito),
                   pct.ribo=mean(pct.ribo))

stats1_rat
plotQC(rat)



# Basic QC filter
rat_filt = doBAsicFiltering(rat)

# correct for contamination with ambient RNA ----


rat_filt_decontx_list = SplitObject(object = rat_filt, split.by = "run10x") %>%   lapply(., function(x) correctAmbientRNA_noEmptyCells(seuratobject = x, ribogenes = ribogenes_rat, mitogenes = mitogenes_rat))

rat_filt_decontx = merge(rat_filt_decontx_list[[1]],y = rat_filt_decontx_list[2:length(rat_filt_decontx_list)], merge.data = TRUE)
rat_filt_decontx

# Identify  Doubletts ----
# refilter as ambient RNA correction might have changed data


rat_filt_decontx2 <- subset(rat_filt_decontx, (pct.mito < param_QC_filter_mito_max) &
                              (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                              (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                              (nCount_RNA < param_QC_filter_RNAcounts_max) &
                              (nCount_RNA > param_QC_filter_RNAcounts_min) ) %>% # & (scrublet_prediction!='doublet')
  doSCTransform(seuratobject = .,experimentIDcol = "run10x", variable.features.n = 4000,   vars.to.regress = "pct.mito") %>%
  CellCycleScoring(s.features = s.genes_rat, g2m.features = g2m.genes_rat,  assay = "SCT") %>%
  RunPCA(features=VariableFeatures(.), assay = "SCT", npcs = 50,verbose=TRUE)%>%
  RunUMAP(dims=1:30,verbose=TRUE,  assay = "SCT" ) %>%
  FindNeighbors(dims=1:30, assay = "SCT", verbose=TRUE)  %>%
  FindClusters( resolution = 0.8) %>%
  FindClusters( resolution = 0.2) %>%
  FindClusters( resolution = 0.4)

rat_filt_decontx2
# plot Cluster ---
plot3clusterings(rat_filt_decontx2)



rat_filt_decontx2 = calcDoubletts(seuratobject = rat_filt_decontx2,
                                  annotation_column = "SCT_snn_res.0.2", # reasonable number of cluster, should roughly match number of cell types
                                  doublets_expected_fromLoading = NULL # corresponding expectation if 8000 cells loaded, see  https://assets.ctfassets.net/an68im79xiti/4tjk4KvXzTWgTs8f3tvUjq/2259891d68c53693e753e1b45e42de2d/CG000183_ChromiumSingleCell3__v3_UG_Rev_C.pdf page 17
)
rat_filt_decontx2

# Filter again including doublets

rat_filt_decontx3 <- subset(rat_filt_decontx2, (pct.mito < param_QC_filter_mito_max) &
                              (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                              (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                              (nCount_RNA < param_QC_filter_RNAcounts_max) &
                              (nCount_RNA > param_QC_filter_RNAcounts_min) &
                              DF.classifications_ConsidHomoDoubl == "Singlet")

rat_filt_decontx3



saveRDS(rat_filt_decontx2, file = here('results/s602_9v2_rat_preprocessed_withDoubletts.RDS'))

# >>HUMAN TRAVAGLINI ----

loaded1 = load(here("data/droplet_normal_lung_seurat_ntiss10x.P2.anno.20191002.RC4.Robj"))
loaded1
ntiss10x.P2.anno
# grep(pattern = "^MT-", rownames(ntiss10x.P2.anno@@assays[["RNA"]]), value = TRUE)
grep("^MT-",dimnames(ntiss10x.P2.anno@raw.data)[[1]],ignore.case = T, value = T)
grep("COX1",dimnames(ntiss10x.P2.anno@raw.data)[[1]],ignore.case = T, value = T)


humanTravaglini = UpdateSeuratObject(ntiss10x.P2.anno)
humanTravaglini

humanTravaglini = DietSeurat(humanTravaglini)


# ensuring basic QC ----
qlist_humanTravaglini0 = venn3(rownames(humanTravaglini), orthologues$`Human gene name`, orthologues$`Human gene stable ID`)
str(qlist_humanTravaglini0)

qlist_humanTravaglini1 = venn3(rownames(humanTravaglini)%>% toupper(), orthologues$`Human gene name`%>% toupper(), orthologues$`Human gene stable ID`%>% toupper())
str(qlist_humanTravaglini1)

grep("^MT-", rownames(humanTravaglini), ignore.case = T, value = T)

mitogenes_humanTravaglini = orthologues[grep('^MT-', `Human gene name`, ignore.case = T), unique(c(`Human gene stable ID`, `Human gene name`)) ] %>% na.omit() %>% intersect(rownames(humanTravaglini) %>% toupper(), .)
mitogenes_humanTravaglini #no mitochondrial genes reported. In the code https://github.com/krasnowlab/HLCA/blob/5dc595048a107692028a9dbd60b5e8a10ea6a660/Annotation%20Patient%202/boilerplate.R only ribosomal percentages are reported, not

ribogenes_humanTravaglini = orthologues[grep('^RP[SL][0-9]*$', `Human gene name`), unique(c(`Human gene stable ID`, `Human gene name`))] %>% na.omit() %>% intersect(rownames(humanTravaglini), .)
ribogenes_humanTravaglini

s.genes_humanTravaglini <- orthologues[`Human gene name` %in% cc.genes$s.genes ,  unique(c(`Human gene stable ID`, `Human gene name`))]  %>% na.omit() %>% intersect(rownames(humanTravaglini), .)
s.genes_humanTravaglini

g2m.genes_humanTravaglini <- orthologues[`Human gene name` %in% cc.genes$g2m.genes, unique(c(`Human gene stable ID`, `Human gene name`))]  %>% na.omit() %>% intersect(rownames(humanTravaglini), .)
g2m.genes_humanTravaglini



orthologues_humanTravaglini = unique(na.omit(c(orthologues$`Human gene stable ID`, orthologues$`Human gene stable ID`)))

# humanTravaglini$pct.mito <- PercentageFeatureSet(humanTravaglini,features = mitogenes_humanTravaglini)
humanTravaglini$pct.ribo =  PercentageFeatureSet(humanTravaglini[rownames(humanTravaglini) %in% orthologues_humanTravaglini,], features = ribogenes_humanTravaglini)

humanTravaglini@meta.data$run10x = paste0("humanTravaglini_",humanTravaglini@meta.data$channel)
table(humanTravaglini@meta.data$run10x)
humanTravaglini@meta.data$species = "humanTravaglini"

# here's some basic statistics.
stats1_humanTravaglini  = humanTravaglini@meta.data %>%
  dplyr::group_by(run10x) %>%
  dplyr::summarise(ncells=dplyr::n(),
                   # ndoublets=sum(scrublet_prediction=='doublet'),
                   nGene=median(nFeature_RNA),
                   nUMI=median(nCount_RNA),
                   # pct.mito=mean(pct.mito),
                   pct.ribo=mean(pct.ribo))

stats1_humanTravaglini
plotQC(humanTravaglini)



# Basic QC filter
humanTravaglini_filt = doBAsicFiltering(humanTravaglini)

# correct for contamination with ambient RNA ----


humanTravaglini_filt_decontx_list = SplitObject(object = humanTravaglini_filt, split.by = "run10x") %>%   lapply(., function(x) correctAmbientRNA_noEmptyCells(seuratobject = x, ribogenes = ribogenes_humanTravaglini, mitogenes = NULL))

humanTravaglini_filt_decontx = merge(humanTravaglini_filt_decontx_list[[1]],y = humanTravaglini_filt_decontx_list[2:length(humanTravaglini_filt_decontx_list)], merge.data = TRUE)

# Identify  Doubletts ----
# refilter as ambient RNA correction might have changed data


humanTravaglini_filt_decontx2 <- subset(humanTravaglini_filt_decontx, (
  nFeature_RNA < param_QC_filter_RNAfeatures_max) &
    # pct.mito < param_QC_filter_mito_max) &
    (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
    (nCount_RNA < param_QC_filter_RNAcounts_max) &
    (nCount_RNA > param_QC_filter_RNAcounts_min) ) %>% # & (scrublet_prediction!='doublet')
  doSCTransform(seuratobject = .,experimentIDcol = "run10x", variable.features.n = 4000,   vars.to.regress = NULL) %>%
  CellCycleScoring(s.features = s.genes_humanTravaglini, g2m.features = g2m.genes_humanTravaglini,  assay = "SCT") %>%
  RunPCA(features=VariableFeatures(.), assay = "SCT", npcs = 50,verbose=TRUE)%>%
  RunUMAP(dims=1:30,verbose=TRUE,  assay = "SCT" ) %>%
  FindNeighbors(dims=1:30, assay = "SCT", verbose=TRUE)  %>%
  FindClusters( resolution = 0.8) %>%
  FindClusters( resolution = 0.2) %>%
  FindClusters( resolution = 0.4)

humanTravaglini_filt_decontx2
# plot Cluster ---
plot3clusterings(humanTravaglini_filt_decontx2)



humanTravaglini_filt_decontx2 = calcDoubletts(seuratobject = humanTravaglini_filt_decontx2,
                                              annotation_column = "SCT_snn_res.0.2", # reasonable number of cluster, should roughly match number of cell types
                                              doublets_expected_fromLoading = NULL # corresponding expectation if 8000 cells loaded, see  https://assets.ctfassets.net/an68im79xiti/4tjk4KvXzTWgTs8f3tvUjq/2259891d68c53693e753e1b45e42de2d/CG000183_ChromiumSingleCell3__v3_UG_Rev_C.pdf page 17
)
humanTravaglini_filt_decontx2

# Filter again including doublets

humanTravaglini_filt_decontx3 <- subset(humanTravaglini_filt_decontx2,
                                        # (pct.mito < param_QC_filter_mito_max) &
                                        (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                                          (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                                          (nCount_RNA < param_QC_filter_RNAcounts_max) &
                                          (nCount_RNA > param_QC_filter_RNAcounts_min) &
                                          DF.classifications_ConsidHomoDoubl == "Singlet")

humanTravaglini_filt_decontx3

table(humanTravaglini_filt_decontx2$run10x)
saveRDS(humanTravaglini_filt_decontx2, file = here('results/s602_9v2_humanTravaglini10x_preprocessed_withDoubletts.RDS'))

# # >>HAMSTER----
hamster = readRDS(here("data/s511_1_diet_hamster.RDS"))
hamster

hamster = DietSeurat(hamster)

table(str_sub(hamster$orig.ident, start = str_length(hamster$orig.ident), end = str_length(hamster$orig.ident)))


hamster = RenameCells(hamster,
                      new.names = paste0(colnames(x = hamster[["RNA"]]),"-", str_sub(hamster$orig.ident, start = str_length(hamster$orig.ident), end = str_length(hamster$orig.ident) ))) # just to make sure unique cell names

# ensuring basic QC ----
qlist_hamster = venn3(rownames(hamster), orthologues$`Golden Hamster gene name`, orthologues$`Golden Hamster gene stable ID`)


mitogenes_hamster = orthologues[grep('^MT-', `Human gene name`), unique(c(`Golden Hamster gene name`,`Golden Hamster gene stable ID`))] %>% na.omit() %>% intersect(rownames(hamster), .)
mitogenes_hamster

ribogenes_hamster = orthologues[grep('^RP[SL][0-9]*$', `Human gene name`), unique(c(`Golden Hamster gene name`,`Golden Hamster gene stable ID`))] %>% na.omit() %>% intersect(rownames(hamster), .)
ribogenes_hamster

s.genes_hamster <- orthologues[`Human gene name` %in% cc.genes$s.genes , unique(c(`Golden Hamster gene name`,`Golden Hamster gene stable ID`))] %>% na.omit() %>% intersect(rownames(hamster), .)
s.genes_hamster

g2m.genes_hamster <- orthologues[`Human gene name` %in% cc.genes$g2m.genes, unique(c(`Golden Hamster gene name`,`Golden Hamster gene stable ID`))] %>% na.omit() %>% intersect(rownames(hamster), .)
g2m.genes_hamster


orthologues_hamster = unique(na.omit(c(orthologues$`Golden Hamster gene name`, orthologues$`Golden Hamster gene stable ID`)))

hamster$pct.mito <- PercentageFeatureSet(hamster[rownames(hamster) %in% orthologues_hamster,],features = mitogenes_hamster)
hamster$pct.ribo =  PercentageFeatureSet(hamster[rownames(hamster) %in% orthologues_hamster,], features = ribogenes_hamster)

hamster@meta.data$run10x = paste0("",hamster@meta.data$orig.ident) # species is as ma present in prefix
table(hamster@meta.data$run10x)
hamster@meta.data$species = "hamster"



# here's some basic statistics.
stats1_hamster  = hamster@meta.data %>%
  dplyr::group_by(run10x) %>%
  dplyr::summarise(ncells=dplyr::n(),
                   # ndoublets=sum(scrublet_prediction=='doublet'),
                   nGene=median(nFeature_RNA),
                   nUMI=median(nCount_RNA),
                   pct.mito=mean(pct.mito),
                   pct.ribo=mean(pct.ribo))

stats1_hamster
plotQC(hamster)



# Basic QC filter
hamster_filt = doBAsicFiltering(hamster)
hh(hamster_filt@assays$RNA@counts)
# Ambient RNA background correction with available emptydroplet data----


todofile_hamster = data.table(run10x = unique(hamster_filt$run10x))
todofile_hamster$raw_h5_10x_fn = c(here("data/GSM4946629_ma-d0-lung-1_raw_feature_bc_matrix.h5"),
                                   here("data/GSM4946630_ma-d0-lung-2_raw_feature_bc_matrix.h5"),
                                   here("data/GSM4946631_ma-d0-lung-3_raw_feature_bc_matrix.h5"))

todofile_hamster$rename_last_2cellidLettersTo = c("-1", "-2", "-3")
todofile_hamster
hamster_filt_decontx_list = SplitObject(object = hamster_filt, split.by = "run10x")
names(hamster_filt_decontx_list)

hamster_filt_decontx_list =  lapply(todofile_hamster$run10x, function(myrun10x) {
  # myrun10x = todofile_hamster$run10x[1]
  myzeile = todofile_hamster[run10x==myrun10x]
  seuratobject = hamster_filt_decontx_list[[myrun10x]]
  seuratobject = correctAmbientRNA_inclEmptyCells(seuratobject=seuratobject,
                                                  raw_h5_10x_fn = myzeile$raw_h5_10x_fn,
                                                  experiment_name <-  paste0("", myrun10x),
                                                  empty_drops_lower = 100,
                                                  ribogenes=ribogenes_hamster,
                                                  mitogenes=mitogenes_hamster ,
                                                  orthologue_genenames = orthologues_hamster,
                                                  rename_last_2cellidLettersTo = myzeile$rename_last_2cellidLettersTo,


                                                  reduce2givenCells=T)

}
)

hamster_filt_decontx = merge(hamster_filt_decontx_list[[1]],y = hamster_filt_decontx_list[2:length(hamster_filt_decontx_list)], merge.data = TRUE)
hamster_filt_decontx


# Clustern, to identifyt  Doubletts ----
# refilter as ambient RNA correction might have changed data


# seuratobject_decontx_withbg = FindVariableFeatures(seuratobject_decontx_withbg, nfeatures = 4000)
hamster_filt_decontx2 <- subset(hamster_filt_decontx, (pct.mito < param_QC_filter_mito_max) &
                                  (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                                  (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                                  (nCount_RNA < param_QC_filter_RNAcounts_max) &
                                  (nCount_RNA > param_QC_filter_RNAcounts_min) ) %>% # & (scrublet_prediction!='doublet')
  doSCTransform(seuratobject = .,experimentIDcol = "run10x", variable.features.n = 4000,   vars.to.regress = "pct.mito")  %>%
  CellCycleScoring(s.features = s.genes_hamster, g2m.features = g2m.genes_hamster,  assay = "SCT")%>%
  RunPCA(features=VariableFeatures(.), assay = "SCT", npcs = 50,verbose=TRUE)%>%
  RunUMAP(dims=1:30,verbose=TRUE,  assay = "SCT" ) %>%
  FindNeighbors(dims=1:30, assay = "SCT", verbose=TRUE)  %>%
  FindClusters( resolution = 0.8) %>%
  FindClusters( resolution = 0.2) %>%
  FindClusters( resolution = 0.4)
hamster_filt_decontx2

# plot Cluster ---

plot3clusterings(hamster_filt_decontx2)


hamster_filt_decontx2 = calcDoubletts(seuratobject = hamster_filt_decontx2,
                                      annotation_column = "SCT_snn_res.0.4", # reasonable number of cluster, should roughly match number of cell types
                                      doublets_expected_fromLoading = NULL # corresponding expectation given number of recovered cells, can be replaced by number of loaded cells, see  https://assets.ctfassets.net/an68im79xiti/4tjk4KvXzTWgTs8f3tvUjq/2259891d68c53693e753e1b45e42de2d/CG000183_ChromiumSingleCell3__v3_UG_Rev_C.pdf
)
hamster_filt_decontx2

# Filter again including doublets


hamster_filt_decontx3 <- subset(hamster_filt_decontx2, (pct.mito < param_QC_filter_mito_max) &
                                  (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                                  (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                                  (nCount_RNA < param_QC_filter_RNAcounts_max) &
                                  (nCount_RNA > param_QC_filter_RNAcounts_min) &
                                  DF.classifications_ConsidHomoDoubl == "Singlet")


hamster_filt_decontx3



saveRDS(hamster_filt_decontx2, file = here('results/s602_9v2_hamster_preprocessed_withDoubletts.RDS'))

# # >>HUMAN CHARITE -----

humanCharite_loaded = readRDS(here("data/s511_1_diet_human_charite.RDS"))
humanCharite_loaded

humanCharite_loaded@meta.data$run10x = paste0("",humanCharite@meta.data$orig.ident) # species is as ma present in prefix
table(humanCharite_loaded@meta.data$run10x)
humanCharite_loaded@meta.data$species = "humanCharite"


humanCharite = DietSeurat(humanCharite_loaded)
#
table(str_sub(humanCharite$orig.ident, start = str_length(humanCharite$orig.ident), end = str_length(humanCharite$orig.ident)))
#
#

# ensuring basic QC ----
qlist_humanCharite = venn3(rownames(humanCharite), orthologues$`Human gene name`, orthologues$`Human gene stable ID`)


mitogenes_humanCharite = orthologues[grep('^MT-', `Human gene name`), unique(c(`Human gene name`,`Human gene stable ID`))] %>% na.omit() %>% intersect(rownames(humanCharite), .)
mitogenes_humanCharite

ribogenes_humanCharite = orthologues[grep('^RP[SL][0-9]*$', `Human gene name`), unique(c(`Human gene name`,`Human gene stable ID`))] %>% na.omit() %>% intersect(rownames(humanCharite), .)
ribogenes_humanCharite

s.genes_humanCharite <- orthologues[`Human gene name` %in% cc.genes$s.genes , unique(c(`Human gene name`,`Human gene stable ID`))] %>% na.omit() %>% intersect(rownames(humanCharite), .)
s.genes_humanCharite

g2m.genes_humanCharite <- orthologues[`Human gene name` %in% cc.genes$g2m.genes, unique(c(`Human gene name`,`Human gene stable ID`))] %>% na.omit() %>% intersect(rownames(humanCharite), .)
g2m.genes_humanCharite


orthologues_humanCharite = unique(na.omit(c(orthologues$`Human gene name`, orthologues$`Human gene stable ID`)))

humanCharite$pct.mito <- PercentageFeatureSet(humanCharite[rownames(humanCharite) %in% orthologues_humanCharite,],features = mitogenes_humanCharite)
humanCharite$pct.ribo =  PercentageFeatureSet(humanCharite[rownames(humanCharite) %in% orthologues_humanCharite,], features = ribogenes_humanCharite)




# here's some basic statistics.
stats1_humanCharite  = humanCharite@meta.data %>%
  dplyr::group_by(run10x) %>%
  dplyr::summarise(ncells=dplyr::n(),
                   # ndoublets=sum(scrublet_prediction=='doublet'),
                   nGene=median(nFeature_RNA),
                   nUMI=median(nCount_RNA),
                   pct.mito=mean(pct.mito),
                   pct.ribo=mean(pct.ribo))

stats1_humanCharite
plotQC(humanCharite)

# Basic QC filter ----
humanCharite_filt = doBAsicFiltering(humanCharite)
hh(humanCharite_filt@assays$RNA@counts)

# Ambient RNA background correction with available emptydroplet data----


todofile_humanCharite = data.table(run10x = unique(humanCharite_filt$run10x))
todofile_humanCharite$raw_h5_10x_fn = c(here("data/lung_700D_control_raw.h5"),
                                        here("data/lung_89C_control_raw.h5"),
                                        here("data/lung_1169Z_control_raw.h5"),
                                        here("data/lung_218V_control_raw.h5"))

todofile_humanCharite$rename_last_2cellidLettersTo = c("-1", "-1", "-1", "-1")
todofile_humanCharite
humanCharite_filt_decontx_list = SplitObject(object = humanCharite_filt, split.by = "run10x")
names(humanCharite_filt_decontx_list)

humanCharite_filt_decontx_list =  lapply(todofile_humanCharite$run10x, function(myrun10x) {
  # humanCharite_filt_decontx_list2 =  lapply(todofile_humanCharite$run10x, function(myrun10x) {

  # myrun10x = todofile_humanCharite$run10x[1]
  myzeile = todofile_humanCharite[run10x==myrun10x]
  seuratobject = humanCharite_filt_decontx_list[[myrun10x]]
  seuratobject = correctAmbientRNA_inclEmptyCells(seuratobject=seuratobject,
                                                  raw_h5_10x_fn = myzeile$raw_h5_10x_fn,
                                                  experiment_name <-  paste0("", myrun10x),
                                                  empty_drops_lower = 100,
                                                  ribogenes=ribogenes_humanCharite,
                                                  mitogenes=mitogenes_humanCharite ,
                                                  orthologue_genenames = orthologues_humanCharite,
                                                  rename_last_2cellidLettersTo = myzeile$rename_last_2cellidLettersTo,
                                                  reduce2givenCells=T,testmode = F)

}
)

# humanCharite_filt_decontx = merge(humanCharite_filt_decontx_list2[[1]],y = humanCharite_filt_decontx_list2[2:length(humanCharite_filt_decontx_list2)], merge.data = TRUE)

humanCharite_filt_decontx = merge(humanCharite_filt_decontx_list[[1]],y = humanCharite_filt_decontx_list[2:length(humanCharite_filt_decontx_list)], merge.data = TRUE)
humanCharite_filt_decontx


# Clustern, to identifyt  Doubletts ----
# refilter as ambient RNA correction might have changed data


# seuratobject_decontx_withbg = FindVariableFeatures(seuratobject_decontx_withbg, nfeatures = 4000)
humanCharite_filt_decontx2 <- subset(humanCharite_filt_decontx, (pct.mito < param_QC_filter_mito_max) &
                                       (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                                       (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                                       (nCount_RNA < param_QC_filter_RNAcounts_max) &
                                       (nCount_RNA > param_QC_filter_RNAcounts_min) ) %>% # & (scrublet_prediction!='doublet')
  doSCTransform(seuratobject = .,experimentIDcol = "run10x", variable.features.n = 4000,   vars.to.regress = "pct.mito") %>%
  CellCycleScoring(s.features = s.genes_humanCharite, g2m.features = g2m.genes_humanCharite,  assay = "SCT")   %>%
  RunPCA(features=VariableFeatures(.), assay = "SCT", npcs = 50,verbose=TRUE)%>%
  RunUMAP(dims=1:30,verbose=TRUE,  assay = "SCT" ) %>%
  FindNeighbors(dims=1:30, assay = "SCT", verbose=TRUE)  %>%
  FindClusters( resolution = 0.8) %>%
  FindClusters( resolution = 0.2) %>%
  FindClusters( resolution = 0.4)
humanCharite_filt_decontx2

# plot Cluster ---

plot3clusterings(humanCharite_filt_decontx2)


humanCharite_filt_decontx2 = calcDoubletts(seuratobject = humanCharite_filt_decontx2,
                                           annotation_column = "SCT_snn_res.0.4", # reasonable number of cluster, should roughly match number of cell types
                                           doublets_expected_fromLoading = NULL # corresponding expectation given number of recovered cells, can be replaced by number of loaded cells, see  https://assets.ctfassets.net/an68im79xiti/4tjk4KvXzTWgTs8f3tvUjq/2259891d68c53693e753e1b45e42de2d/CG000183_ChromiumSingleCell3__v3_UG_Rev_C.pdf
)
humanCharite_filt_decontx2

# Filter again including doublets


humanCharite_filt_decontx3 <- subset(humanCharite_filt_decontx2, (pct.mito < param_QC_filter_mito_max) &
                                       (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                                       (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                                       (nCount_RNA < param_QC_filter_RNAcounts_max) &
                                       (nCount_RNA > param_QC_filter_RNAcounts_min) &
                                       DF.classifications_ConsidHomoDoubl == "Singlet")


humanCharite_filt_decontx3

saveRDS(humanCharite_filt_decontx2, file = here('results/s602_9v2_humanCharite_preprocessed_withDoubletts.RDS'))




# # >>HUMAN SS2 -----
# drop seq sequencing from Travaligni

loaded1 = load(here("data/facs_normal_lung_blood_seurat_ntiss.P2.anno.gencode.20200616.RC4.Robj"))
loaded1
ntiss.P2.anno.gencode

# grep(pattern = "^MT-", rownames(ntiss.P2.anno.gencode@@assays[["RNA"]]), value = TRUE)
grep("^MT-",dimnames(ntiss.P2.anno.gencode@raw.data)[[1]],ignore.case = T, value = T)
grep("COX1",dimnames(ntiss.P2.anno.gencode@raw.data)[[1]],ignore.case = T, value = T)


humanSS2 = UpdateSeuratObject(ntiss.P2.anno.gencode)
humanSS2

feats <- c("nFeature_RNA", "nCount_RNA")
VlnPlot(humanSS2, group.by = "orig.ident", features = feats, pt.size = 0.1)
summary(humanSS2$nReads)

humanSS2 = DietSeurat(humanSS2)


# ensuring basic QC ----
qlist_humanSS20 = venn3(rownames(humanSS2), orthologues$`Human gene name`, orthologues$`Human gene stable ID`)
str(qlist_humanSS20)

qlist_humanSS21 = venn3(rownames(humanSS2)%>% toupper(), orthologues$`Human gene name`%>% toupper(), orthologues$`Human gene stable ID`%>% toupper())
str(qlist_humanSS21)

grep("^MT-", rownames(humanSS2), ignore.case = T, value = T)

mitogenes_humanSS2 = orthologues[grep('^MT-', `Human gene name`, ignore.case = T), unique(c(`Human gene stable ID`, `Human gene name`)) ] %>% na.omit() %>% intersect(rownames(humanSS2) %>% toupper(), .)
mitogenes_humanSS2 # mitochondrial genes reported, in contrast to the 10x data

ribogenes_humanSS2 = orthologues[grep('^RP[SL][0-9]*$', `Human gene name`), unique(c(`Human gene stable ID`, `Human gene name`))] %>% na.omit() %>% intersect(rownames(humanSS2), .)
ribogenes_humanSS2

s.genes_humanSS2 <- orthologues[`Human gene name` %in% cc.genes$s.genes ,  unique(c(`Human gene stable ID`, `Human gene name`))]  %>% na.omit() %>% intersect(rownames(humanSS2), .)
s.genes_humanSS2

g2m.genes_humanSS2 <- orthologues[`Human gene name` %in% cc.genes$g2m.genes, unique(c(`Human gene stable ID`, `Human gene name`))]  %>% na.omit() %>% intersect(rownames(humanSS2), .)
g2m.genes_humanSS2



orthologues_humanSS2 = unique(na.omit(c(orthologues$`Human gene stable ID`, orthologues$`Human gene stable ID`)))

humanSS2$pct.mito <- PercentageFeatureSet(humanSS2,features = mitogenes_humanSS2)

humanSS2$pct.ribo =  PercentageFeatureSet(humanSS2[rownames(humanSS2) %in% orthologues_humanSS2,], features = ribogenes_humanSS2)

humanSS2@meta.data$run10x = paste0("humanSS2_",humanSS2@meta.data$orig.ident)
table(humanSS2@meta.data$run10x)
humanSS2@meta.data$species = "humanSS2"

# here's some basic statistics.
stats1_humanSS2  = humanSS2@meta.data %>%
  dplyr::group_by(run10x) %>%
  dplyr::summarise(ncells=dplyr::n(),
                   # ndoublets=sum(scrublet_prediction=='doublet'),
                   nGene=median(nFeature_RNA),
                   nUMI=median(nCount_RNA),
                   pct.mito=mean(pct.mito),
                   pct.ribo=mean(pct.ribo))

stats1_humanSS2
plotQC(humanSS2)



# Basic QC filter
humanSS2_filt = doBAsicFiltering(humanSS2, QC_filter_RNAcounts_max = 4000000, QC_filter_RNAcounts_min = 50000)

# correct for contamination with ambient RNA ----
# not required

# Identify  Doubletts ----
# not required, but calculated as sensitivity analysis


humanSS2_filt <- humanSS2_filt %>% # & (scrublet_prediction!='doublet')
  doSCTransform(seuratobject = .,experimentIDcol = "run10x", variable.features.n = 4000,   vars.to.regress = NULL) %>%
  CellCycleScoring(s.features = s.genes_humanSS2, g2m.features = g2m.genes_humanSS2,  assay = "SCT") %>%
  RunPCA(features=VariableFeatures(.), assay = "SCT", npcs = 50,verbose=TRUE)%>%
  RunUMAP(dims=1:30,verbose=TRUE,  assay = "SCT" ) %>%
  FindNeighbors(dims=1:30, assay = "SCT", verbose=TRUE)  %>%
  FindClusters( resolution = 0.8) %>%
  FindClusters( resolution = 0.2) %>%
  FindClusters( resolution = 0.4)

humanSS2_filt
# plot Cluster ---
plot3clusterings(humanSS2_filt)



humanSS2_filt_decontx2 = calcDoubletts(seuratobject = humanSS2_filt,
                                       annotation_column = "SCT_snn_res.0.4", # reasonable number of cluster, should roughly match number of cell types
                                       doublets_expected_fromLoading = 0.005 # in SS2 ,I expect little amount of doublets
)
humanSS2_filt_decontx2
table(humanSS2_filt_decontx2$DF.classifications_ConsidHomoDoubl)

DimPlot(humanSS2_filt_decontx2, group.by = "free_annotation") %>% plotly::ggplotly()
# Filter again including doublets

humanSS2_filt_decontx3 <- subset(humanSS2_filt_decontx2,
                                 (pct.mito < param_QC_filter_mito_max) &
                                   (nFeature_RNA < param_QC_filter_RNAfeatures_max) &
                                   (nFeature_RNA > param_QC_filter_RNAfeatures_min) &
                                   (nCount_RNA < 4000000) &
                                   (nCount_RNA > param_QC_filter_RNAcounts_min) &
                                   DF.classifications_ConsidHomoDoubl == "Singlet")

humanSS2_filt_decontx3


saveRDS(humanSS2_filt_decontx2, file = here('results/s602_9v2_humanTravagliniSS2_preprocessed_withDoubletts.RDS'))

showNA(orthologues, showAllNoNA = F)
fwrite(orthologues, here("results/s602_9v2_orthologues.txt.gz"))


finalizeSkript()


#
# require(toolboxH)
# require(here)
# require(rmarkdown)
#
# render(here("scripts/s602_9v2_preprocess_singlestudies.R"),encoding="UTF-8")
