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

packages2load = c("magrittr","knitr", "ggplot2", "scales", "ggthemes", "assertr","Seurat", 'paletteer',"pheatmap", "cowplot", "here", "plotly", "ggrepel", "ggforce", "dplyr")

new.packages <- packages2load[!(packages2load %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

for(i in packages2load) {

  suppressPackageStartupMessages(library(i, character.only = TRUE))

}


# CREDITS to https://stackoverflow.com/questions/48184645/how-can-i-put-the-labels-outside-of-piechart

# # LOAD DATA ----
#
lunge_merged = readRDS(here("results/s1_3_lunge_merged.RDS")) # processed integrated data, available at https://www.health-atlas.de/studies/54


# # CREATE DATA ----

lunge_merged.meta.data = data.table(lunge_merged@meta.data, keep.rownames = T)

tissue2groupm_counted = lunge_merged.meta.data[,.N, .(cluster_seurat_v2,  species)]


# tissue2groupm_counted[, cluster_seurat_v2:= factor(cluster_seurat_v2, levels = reihenfolge)]

tissue2groupm_counted[,prozi := N/sum(N), species]
tissue2groupm_counted



setorder(tissue2groupm_counted, species,cluster_seurat_v2)
tissue2groupm_counted2 <- tissue2groupm_counted %>% group_by(species) %>%
  mutate(end = 2 * pi * cumsum(N)/sum(N),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))

tissue2groupm_counted2 = data.table(tissue2groupm_counted2)


levels(tissue2groupm_counted2$cluster_seurat_v2)

# tissue2groupm_counted2[,species:= factor(species, levels = c('Hamster', "Mouse" , "HumanK", "HumanT"))]

tissue2groupm_counted2[,cluster_seurat_v22 :=  str_replace_all(cluster_seurat_v2, "Proliferating", "Prolif.") %>% str_wrap(., width = 16)]
setorder(tissue2groupm_counted, cluster_seurat_v2)
tissue2groupm_counted2[,cluster_seurat_v22 := factor(cluster_seurat_v22, levels = unique(cluster_seurat_v22))]



# # PLOT ----
nogpalette3 = c("#CB769E", "#DE639A", "#A85C85", "#7C6A0A", "#4F6D7A", "#368F8B",
                "#246A73", "#5CC1BC", "#62C370", "#F7C548", "#F97E44", "#FB3640",
                "#B7245C", "#0D3B66", "#3E2F5B", "#B2675E", "#644536")


pieplot = ggplot(tissue2groupm_counted2,aes( fill = cluster_seurat_v22, color = cluster_seurat_v22)) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start, end = end)
               # , col = "white"
  )+
  facet_wrap(~species, nrow = 1) +
  coord_fixed()  +
  # geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = cluster_seurat_v2,
  # hjust = hjust, vjust = vjust),size=2)+ #,min.segment.length=0.1, force_pull = 0.1000) +

  scale_x_continuous(limits = c(-1.01, 1.01),  # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1.01, 1.01),      # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL) +
  # guides(fill = "none", color = "none") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE, title = element_blank()),
         col=guide_legend(nrow=2,byrow=TRUE, title = element_blank())) +
  theme_pander(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 13.5),
        legend.key.size = unit(1,"cm"),
        panel.spacing = unit(1.2, "cm"),
        strip.text.x =  element_text(size = 18)) +
  scale_color_manual(values =nogpalette3 ) +
  scale_fill_manual(values = nogpalette3)

pieplot

jpeg(here("results/s2_1_figure1B_plot_celltype_proportions.jpeg"),width = 18,height = 5,units = "in", quality = 100, res = 150)
print(pieplot)
dev.off()

# # FINALIZE SCRIPT----

finalizeSkript()
