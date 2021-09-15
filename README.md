
# A pulmonologist’s guide to perform and analyse cross-species single-lung-cell transcriptomics 

Single-cell ribonucleic acid sequencing (scRNA-seq) is becoming widely employed to study biological processes in great depth. Its name describes its key advantage, the resolution of transcriptomes at single-cell level.
In the Publication mentioned above, we exemplarily demonstrate integration and qualitative assessment of scRNA-Seq-data across three species: human (Homo sapiens), hamster (Mesocricetus auratus), and mice (Mus musculus). We used quality-filtered scRNA-seq data of lungs from 4 different sources (Figure 1) using up to date analysis workflows.
All cross-species analysis code of this example including expression of canonical celltype-specific genes is made available via this Github page. For the creation of Figure 1 further processing and alignment was performed with Inkscape (Inkscape Community; https://inkscape.org/de/; gitlab.com/inkscape/inkscape). 

Input files as well as the processed, integrated result file are available as Seurat-R-object at the [Leipzig Health Atlas](https://www.health-atlas.de/studies/54). 



![Figure Method 210830_v5 3](https://user-images.githubusercontent.com/73164857/132707879-fd342a69-1ba5-4e44-bc23-9d0d9f737a22.PNG)
FIGURE 1: 
Cross-species comparison of sequenced lung cells. A) Schematic single-cell ribonucleic acid sequencing workflow, created with Biorender.com. B) Pie charts of lung cell frequencies in indicated species estimated from scRNA-Seq data. C) Uniform manifold approximation and projection (UMAP) plot of identified cell populations across species. D) Feature plots indicating the stress response of single lung cells across species based on the 12 stress related genes _ATF3,  BTG2, DUSP1, EGR1, FOS, FOSB, HSPA8, HSPB1, IER3, JUN, JUNB,_ and _JUND_. Prolif.: proliferating, AT1: alveolar epithelial cells type 1, AT2: alveolar epithelial cells type 2, ly Endothelial: lymphatic endothelial cells. Samples: Mouse, 3205 cells of two whole lungs pooled prior analysis. Hamster, 10564 cells from lobus caudalis of three hamsters from [Nouailles-Kursar et al.]( https://doi.org/10.1038/s41467-021-25030-7). Human Charité, 18487 cells from four fresh lung explants. [Human Travaglini et al.](https://doi.org/10.1038/s41586-020-2922-4), 28793 cells of patient 2.

# Supplementary Figures to the manuscript 

![Dotplot_Marker_All in One](https://github.com/GenStatLeipzig/pulmonologists_interspecies_scRNA/blob/main/results/s5_1_canonical_markers.jpeg)


FIGURE S1: 
Dotplot of cross-species comparison of RNA marker genes. The point size indicates the percentage of cells in a certain population expressing the indicated marker gene. The colour indicates the dataset, and the transparancy indicates the expression level. Hamster: red, Mouse: blue, Human Charité: black, Human Travaglini et al.: green. Average expression is indicated by colour intensity. 

# Supplementary Informations on sample material

Hamster: Dataset consists of lung cells isolated from lobus caudalis of two female and one male Syrian hamsters (Mesocricetus auratus; breed RjHan:AURA, Janvier Labs, France) at 10 to 12 weeks of age. It includes expression information for 35604 expression features for a total of 10564 cells resulting from three experimental batches. All experiments involving animals were approved by institutional and governmental authorities (Freie Universität Berlin and LaGeSo Landesamt für Gesundheit und Soziales Berlin, Germany). For further Information see: https://doi.org/10.1038/s41467-021-25030-7

Mouse: Dataset consists of pulmonary cells isolated from whole lung tissue. The single cell suspension of two 8 to 10 weeks old C57BL/6J mice was pooled prior performing single cell barcoding and library construction according to manufacturer’s instructions (10X Genomics). All experiments involving animals were approved by institutional and governmental authorities (Charité Universitätsmedizin Berlin and LaGeSo Landesamt für Gesundheit und Soziales Berlin, Germany). Alignment and quantification was done using cellranger against the human GRCh38 genome. Low-quality cells were removed by filtering out cells not having less than 1000 and not more than 35000 RNA counts quantified as unique molecular identifier (UMIs), not having at least 300 different genes expressed and not having more than 30% mitochondrial genes expressed. Data was normalized applying standard Seurat 4.0.1 SCTransform workflow, the mouse data was available as a single batch. The dataset includes expression information for 45046 expression features for a total of 3205 cells 
.

Human Charité: Dataset consists of lung cells isolated from fresh lung explants of four 65 - 85 years old patients (3 male, 1 female), suffering from lung carcinoma who underwent lung resection. Only tumor-free peripheral lung tissue was used. The study was approved by the ethics committee at the Charité Universitätsmedizin Berlin (projects EA2/079/13). Written informed consent was obtained from all patients. Single cell barcoding and library construction was done according to the manufacturer’s instructions (10X Genomics, run in four experimental batches). Alignment and quantification was done using cellranger against the human GRCh38 genome. Standard Seurat workflow was performed, excluding cells with 10% or more mitochondrial genes and less than 500 genes, batch effect with were removed using Seurat::IntegrateData(), and doublets were identified using [scrublet](https://github.com/swolock/scrublet). The dataset includes expression information for 44621 expression features for a total of 18487 cells.

Human Travaglini et al.: Sequenced lung cell data of patient 2 published at https://www.synapse.org by Travaglini et al. (https://doi.org/10.1038/s41586-020-2922-4) was used. According to the original publication the cells were isolated from freshly resected lung tissue of a 46-year-old male. The dataset used here includes expression information for 26485 expression features for a total of 28793 cells originating from eight experimental batches.
