rm(list = ls())
require(toolboxH)
require(here)
require(Seurat)
require(ggplot2)
require(ggthemes)
require(dplyr)
require(tidyr)
require(patchwork)


# >>ratOri ---- 
ratOri = readRDS(file = here("data/s101_2_rat_DS024_DS025.RDS"))
ratOri

orthologues = fread(here("results/s602_9v2_orthologues.txt.gz"))

param_QC_filter_mito_max 	=30 # ABX 30, human 10
param_QC_filter_RNAfeatures_min 	=300
param_QC_filter_RNAfeatures_max 	= 6000
param_QC_filter_RNAcounts_min 	=1000
param_QC_filter_RNAcounts_max	=35000
# replace "" with NA in Gene names----
showNA(orthologues, showAllNoNA = F)
for(i in c( "Rat gene name", "Rat gene stable ID")) {
  orthologues[get(i)=="",(i):=NA]
}
showNA(orthologues, showAllNoNA = F)


# ensuring basic QC ----
qlist_rat0 = venn3(rownames(ratOri), orthologues$`Rat gene name`, orthologues$`Rat gene stable ID`)
str(qlist_rat0)

qlist_rat1 = venn3(rownames(ratOri)%>% toupper(), orthologues$`Rat gene name`%>% toupper(), orthologues$`Rat gene stable ID`%>% toupper())
str(qlist_rat1)
table(str_sub(orthologues$`Rat gene stable ID`, 1,6))

# change gene names to upper case to overcome discrepancies with orthologue names
orthologues$`Rat gene name uc` = toupper(orthologues$`Rat gene name` )

ratOri = Seurat.utils::RenameGenesSeurat(obj = ratOri, newnames = toupper(rownames(ratOri))) # # from https://github.com/vertesy/Seurat.utils/
ratOri

qlist_rat02 = venn3(rownames(ratOri), orthologues$`Rat gene name uc`, orthologues$`Rat gene stable ID`)
str(qlist_rat02)

mitogenes_rat = orthologues[grep('^MT-', `Human gene name`), unique(c(`Rat gene stable ID`, `Rat gene name uc`)) %>% toupper()] %>% na.omit() %>% intersect(rownames(ratOri)%>% toupper(), .)
mitogenes_rat = rownames(ratOri)[rownames(ratOri)%>% toupper() %in% mitogenes_rat]
mitogenes_rat

orthologues_rat = unique(na.omit(c(orthologues$`Rat gene name uc`, orthologues$`Rat gene stable ID`)))

mitogenes_rat %in% rownames(ratOri[rownames(ratOri) %in% orthologues_rat,])

ratOri$pct.mito <- PercentageFeatureSet(ratOri[rownames(ratOri) %in% orthologues_rat,],features = mitogenes_rat)

ratOri@meta.data$run10x = paste0("rat_",ratOri@meta.data$orig.ident)
table(ratOri@meta.data$run10x)
ratOri@meta.data$species = "rat"


# plot for the methods summary figure ----
param_QC_filter_mito_max 	=30 # ABX 30, human 10


qccols = c('run10x', 'nCount_RNA', 'nFeature_RNA',  'pct.mito')
plotdata = ratOri@meta.data %>% dplyr::select(qccols[qccols %in% names(ratOri@meta.data)]) %>% as.data.table()

plotdatam =  melt(plotdata, id.vars = 'run10x', variable.name = 'metric')

plotdatam[, metric2 := ifelse(metric == 'nCount_RNA', "Counts\nper cell", 
                              ifelse(metric == "nFeature_RNA", "Transcripts\nper cell", 
                                     ifelse(metric== "pct.mito", "% Mitochondrial\ngenes per cell", metric)))]


plotdatam[, run10x2 := factor(run10x) %>% as.numeric %>% paste0("Batch ",.)]
plotdatam[, .N, .(run10x,run10x2)]
plotdatam[, .N, .(metric, metric2)]

linesdat = data.table(metric2 = unique(plotdatam$metric2)[c(1,1,2,2,3)])

linesdat[,hline := c(param_QC_filter_RNAcounts_min, 
                     param_QC_filter_RNAcounts_max, 
                     param_QC_filter_RNAfeatures_min, 
                     param_QC_filter_RNAfeatures_max,
                     param_QC_filter_mito_max)]

p1 = ggplot(plotdatam, 
            aes(x=run10x2,y=value,fill=run10x2)) + 
  geom_violin() + 
  facet_wrap(~metric2,ncol=4,scales='free_y') + 
  theme_minimal(base_size = 16)+
  xlab('Two batches rat lungs, public data')+
  theme(
    axis.text.x = element_blank()
    # axis.title.x=element_blank(),
    # axis.text.x = element_text(angle = 90, vjust = 0.4, hjust  =0)
  ) +
  
  guides(fill  = "none") +
  geom_hline(data = linesdat,aes(yintercept = hline), col = "orange2", lwd = 2)


p1

jpeg(here("results/s608a_1_methodsfigure_prepro.jpeg"), 6, 2.5, units = "in", quality = 100, res = 150)
p1
dev.off()


finalizeSkript()
