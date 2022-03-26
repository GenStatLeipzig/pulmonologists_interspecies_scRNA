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
set.seed(2712)
randomgenes1k = sample(rownames(lunge_merged_processed), 500)
str(randomgenes1k)
matrix_RNA = lunge_merged_processed@assays$RNA[randomgenes1k,]
matrix_RNAnum_pre = as.numeric(matrix_RNA)
matrix_RNAnum = matrix_RNAnum_pre[matrix_RNAnum_pre>0]

hist(matrix_RNAnum, breaks = 200)
hist(matrix_RNAnum, breaks = 200, ylim = c(0, 100000))
hist(matrix_RNAnum[ matrix_RNAnum<max(matrix_RNAnum)/5 ], breaks = 200)


matrix_SCT = lunge_merged_processed@assays$SCT[randomgenes1k,]
matrix_SCTnum_pre = as.numeric(matrix_SCT)
matrix_SCTnum = matrix_SCTnum_pre[matrix_SCTnum_pre>0]
hist(matrix_SCTnum, breaks = 200)
hist(matrix_SCTnum[ matrix_SCTnum<max(matrix_SCTnum)/5 ], breaks = 200)

require("ggforce")

ggplot(iris, aes(Petal.Length, Petal.Width, colour = Species)) +
  geom_point() +
  facet_zoom(xlim = c(2, 4))

require(scales)
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
library(scales)



matrix_RNAnum_dt = data.table(matrix_RNAnum)
p1 = ggplot(matrix_RNAnum_dt, aes(matrix_RNAnum)) + geom_histogram(bins = 400) + xlab("RNA counts")+ylab("Observations\n(500 random genes)")+
  facet_zoom(xlim = c(1, 100), horizontal = F)+ theme_solarized()+scale_y_continuous(label=label_comma())

p1

matrix_SCTnum_dt = data.table(matrix_SCTnum)
p0 = ggplot(matrix_SCTnum_dt, aes(matrix_SCTnum)) + geom_histogram(bins = 200) + xlab("normalized and\n transformed counts")+
  facet_zoom(xlim = c(0.5, 1.5), horizontal = F) + theme_solarized()+scale_y_continuous(label=label_comma()) +ylab("Observations\n(500 random genes)")
  
p0



jpeg2(here("results/s607_3_no_transform.jpeg"), 3,3)
p1
dev.off()

jpeg2(here("results/s607_3_SCT_transformed.jpeg"), 3,3)
p0
dev.off()

# nochmal

set.seed(2712)
randomgenes1k = sample(rownames(lunge_merged_processed), 100)
str(randomgenes1k)
matrix_RNA = lunge_merged_processed@assays$RNA[randomgenes1k,]
matrix_RNAnum_pre = as.numeric(matrix_RNA)
matrix_RNAnum = matrix_RNAnum_pre[matrix_RNAnum_pre>0]

hist(matrix_RNAnum, breaks = 200)
hist(matrix_RNAnum, breaks = 200, ylim = c(0, 100000))
hist(matrix_RNAnum[ matrix_RNAnum<max(matrix_RNAnum)/5 ], breaks = 200)


matrix_SCT = lunge_merged_processed@assays$SCT[randomgenes1k,]
matrix_SCTnum_pre = as.numeric(matrix_SCT)
matrix_SCTnum = matrix_SCTnum_pre[matrix_SCTnum_pre>0]
hist(matrix_SCTnum, breaks = 200)
hist(matrix_SCTnum[ matrix_SCTnum<max(matrix_SCTnum)/5 ], breaks = 200)

library(scales)

matrix_RNAnum_dt = data.table(matrix_RNAnum)
p1v2 = ggplot(matrix_RNAnum_dt, aes(matrix_RNAnum)) + geom_histogram( binwidth=1) + xlab("RNA counts")+ylab("Observations\n(100 random genes)")+
theme_minimal_grid()+scale_y_continuous(label=label_comma())

p1v2

matrix_SCTnum_dt = data.table(matrix_SCTnum)
p0v2 = ggplot(matrix_SCTnum_dt, aes(matrix_SCTnum)) + geom_histogram( binwidth = 0.5) + xlab("normalized and\n transformed counts")+
  theme_minimal_grid()+scale_y_continuous(label=label_comma()) +ylab("Observations\n(100 random genes)") + xlim(c(0, 5))

p0v2



jpeg2(here("results/s607_3_no_transform.jpeg"), 3,3)
p1v2
dev.off()

jpeg2(here("results/s607_3_SCT_transformed.jpeg"), 3,3)
p0v2
dev.off()
