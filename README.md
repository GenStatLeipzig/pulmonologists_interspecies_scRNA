
# A pulmonologist’s guide to perform and analyse cross-species single-lung-cell transcriptomics 

Single-cell ribonucleic acid sequencing (scRNA-seq) is becoming widely employed to study biological processes in great depth. Its name describes its key advantage, the resolution of transcriptomes at single-cell level.
In the Publication mentioned above, we exemplarily demonstrate integration and qualitative assessment of scRNA-Seq-data across three species: human (Homo sapiens), hamster (Mesocricetus auratus), and mice (Mus musculus). We used quality-filtered scRNA-seq data of lungs from 4 different sources (Figure 1) using up to date analysis workflows.
All cross-species analysis code of this example including expression of canonical celltype-specific genes is made available via this Github page. For the creation of Figure 1 further processing and alignment was performed with Inkscape (Inkscape Community; https://inkscape.org/de/; gitlab.com/inkscape/inkscape)



![Figure Method 210830_v5 3](https://user-images.githubusercontent.com/73164857/132707879-fd342a69-1ba5-4e44-bc23-9d0d9f737a22.PNG)
FIGURE 1: 
Cross-species comparison of sequenced lung cells. A) Schematic single-cell ribonucleic acid sequencing workflow, created with Biorender.com. B) Pie charts of lung cell frequencies in indicated species. C) Uniform manifold approximation and projection (UMAP) plot of identified cell populations across species. D) Feature plots indicating the stress response of single lung cells across species based on the 14 stress related genes ATF3, BTG1, BTG2, DUSP1, EGR1, FOS, FOSB, HSP90AB1, HSPA8, HSPB1, IER3, JUN, JUNB, and JUND. Prolif.: proliferating, AT1: alveolar epithelial cells type 1, AT2: alveolar epithelial cells type 2, ly Endothelial: lymphatic endothelial cells.




# Supplementary Figures to the manuscript 

![Dotplot_Marker_All in One](https://github.com/GenStatLeipzig/pulmonologists_interspecies_scRNA/blob/main/results/s5_1_canonical_markers.jpeg)


FIGURE S1: 
Dotplot of cross-species comparison of RNA marker genes. The point size indicates the percentage of cells in a certain population expressing the indicated marker gene. The colour indicates the dataset. Hamster: red, Mouse: blue, Human Charité: green, Human Travaglini et al.: black. Average expression is indicated by colour intensity. 

# Supplementary Informations on sample material

Hamster: Dataset consists of lung cells isolated from lobus caudalis of two female and one male Syrian hamsters (Mesocricetus auratus; breed RjHan:AURA, Janvier Labs, France) at 10 to 12 weeks of age. All experiments involving animals were approved by institutional and governmental authorities (Freie Universität Berlin and LaGeSo Landesamt für Gesundheit und Soziales Berlin, Germany). For further Information see: https://doi.org/10.1038/s41467-021-25030-7

Mouse: Dataset consists of pulmonary cells isolated from whole lung tissue. The single cell suspension of two 8 to 10 weeks old C57BL/6J mice was pooled prior performing single cell barcoding and library construction according to manufacturer’s instructions (10X Genomics). All experiments involving animals were approved by institutional and governmental authorities (Charité Universitätsmedizin Berlin and LaGeSo Landesamt für Gesundheit und Soziales Berlin, Germany).

Human Charité: Dataset consists of lung cells isolated from fresh lung explants of four patients suffering from lung carcinoma who underwent lung resection. Only tumor-free peripheral lung tissue was used. The study was approved by the ethics committee at the Charité Universitätsmedizin Berlin. Written informed consent was obtained from all patients.

Human Travaglini et al.: Sequenced lung cell data of patient 2 published at https://www.synapse.org by Travaglini et al. (https://doi.org/10.1038/s41586-020-2922-4) was used. According to the original publication the cells were isolated from freshly resected lung tissue of a 46-year-old male.
