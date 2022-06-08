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

jpeg2 = function(..., res = 150) jpeg(..., units = "in", quality = 100, res = res)

# # LOADING----
lunge_merged_processed = readRDS( here("results/s603_2_seurat_lunge_merged_21-12-11.RDS"))

p0 <- DimPlot(lunge_merged_processed, reduction = "umap", group.by = "species", label = T) + NoLegend() + ggtitle("") + xlab("UMAP 1")  + ylab("UMAP 2")
p0
p0data = p0$data %>% data.table(., keep.rownames=T)
p0data

labeltext = p0data[, .(UMAP_1 = -10,
                       UMAP_2 = seq(5, 13, (13-5)/5),
                   species=unique(species))
                   ]
labeltext
require(cowplot)
p0v2 = ggplot(p0data, aes(UMAP_1, UMAP_2, col = species, label = species)) + geom_point(size = 0.3, alpha = 0.2) + geom_text(data = labeltext,fontface = "bold", alpha = 0.8 ) + theme_cowplot() + guides(col = "none", label = "none")+ ggtitle("") + xlab("UMAP 1")  + ylab("UMAP 2")
p0v2

pdf(here("results/s607_3_umap_with_integration.pdf"), 3.5,3.5)
p0v2
dev.off()

jpeg2(here("results/s607_3_umap_with_integration.jpeg"), 3.5,3.5)
p0v2
dev.off()



# Without integration

lunge_merged = DietSeurat(lunge_merged_processed)

lunge_merged@active.assay
lunge_merged@active.assay = "RNA"
lunge_merged@active.assay


lunge_merged2 <- lunge_merged %>%
  NormalizeData(normalization.method='LogNormalize',scale.factor=10000,verbose=TRUE) %>%
  FindVariableFeatures(selection.method='vst',nfeatures=2000,verbose=TRUE) %>%
  ScaleData(vars.to.regress=c('nCount_RNA'),verbose=TRUE) %>%
  # CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) %>%
  RunPCA(features=VariableFeatures(.),verbose=TRUE) %>%
  FindNeighbors(dims=1:50,verbose=TRUE) %>%
  FindClusters(resolution = 0.3) %>%
  FindClusters(resolution = 0.4) %>%
  RunUMAP(dims=1:50,verbose=TRUE)

p1 <- DimPlot(lunge_merged2, reduction = "umap", group.by = "species", label = T) + NoLegend() + ggtitle("") + xlab("UMAP 1")  + ylab("UMAP 2")
p1
p1data = p1$data %>% data.table(., keep.rownames=T)
p1data

labeltext = p1data[, .(UMAP_1 = -12.5,
                       UMAP_2 = seq(9, 17, (17-9)/5),
                       species=unique(species))
]
labeltext
require(cowplot)
p1v2 = ggplot(p1data, aes(UMAP_1, UMAP_2, col = species, label = species)) + geom_point(size = 0.3, alpha = 0.2) + geom_text(data = labeltext,fontface = "bold" , alpha = 0.8) + theme_cowplot() + guides(col = "none", label = "none")+ ggtitle("") + xlab("UMAP 1")  + ylab("UMAP 2")
p1v2

pdf(here("results/s607_3_umap_NO_integration.pdf"), 3.5,3.5)
p1v2
dev.off()

jpeg2(here("results/s607_3_umap_NO_integration.jpeg"), 3.5,3.5)
p1v2
dev.off()


## variance stabilistation ----

set.seed(1802)
randomgenes1k = sample(rownames(lunge_merged_processed), 20)
str(randomgenes1k)

matrix_RNA = lunge_merged_processed@assays$RNA[randomgenes1k,]
dim(matrix_RNA)
matrix_RNAnum_pre = as.data.table(matrix_RNA %>% as.matrix %>% t)
matrix_RNAnum = melt(matrix_RNAnum_pre)
ngen=10
p1v3 = ggplot(matrix_RNAnum[variable %in% unique(variable)[1:ngen]&value>0], aes(variable, value, fill = variable)) + geom_boxplot(fill='red', col ="red") + xlab(paste0(ngen, "randomly\nchosen genes"))+ylab("Raw\nCounts per cell")+
  theme_minimal_grid()+scale_y_continuous(label=label_comma())+
  theme(
    # axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.4)+
    axis.text.x = element_blank()
  )


p1v3

matrix_SCT = lunge_merged_processed@assays$SCT[randomgenes1k,]
dim(matrix_SCT)
matrix_SCTnum_pre = as.data.table(matrix_SCT %>% as.matrix %>% t)
matrix_SCTnum = melt(matrix_SCTnum_pre)
matrix_SCTnum = matrix_SCTnum[ value>0]
matrix_SCTnum[, value:= scale(value ), variable]
p0v3 = ggplot(matrix_SCTnum[variable %in% unique(variable)[1:ngen]], aes(variable, value, fill = variable)) + geom_boxplot(fill='red', col ="red") + xlab(paste0(ngen, " randomly\nchosen genes"))+ylab("SCT-transformed\nscaled Counts per cell")+
  theme_minimal_grid()+scale_y_continuous(label=label_comma())+
  theme(
    # axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.4)+
    axis.text.x = element_blank()
  )
p0v3

p1v3 + p0v3


jpeg2(here("results/s607_3v3_no_transform.jpeg"), 3.3,3)
p1v3
dev.off()

jpeg2(here("results/s607_3v3_SCT_transformed.jpeg"), 3.3,3)
p0v3
dev.off()
