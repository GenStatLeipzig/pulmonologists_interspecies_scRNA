rm(list = ls())

library(SingleCellSignalR)
library(Seurat)
library(generics)
library(here)
require(toolboxH)
require(stringr)
require(ggplot2)
require(ggthemes)
require(scales)
require(cowplot)

plotSankey = function (comparematrix, spalte4color, gap.width = 0.5, colorlines = T, 
          alpha = 0.65, showblocks = F, x_axis_size = 0.1, mymargins = c(2, 1, 1, 1),...) 
{
  comparematrix = data.table::as.data.table(data.table::copy(comparematrix))
  matr_names = names(comparematrix)
  matr_count = comparematrix[, .N, mget(matr_names)]
  for (i in matr_names) {
    matr_count[is.na(get(i)), `:=`((i), "NA")]
  }
  matr_count[, `:=`(orinum, .I)]
  matr_count[, `:=`(col1_factor_level, as.numeric(as.factor(get(matr_names[1]))) %>% 
                      formatC(., width = 6, format = "d", flag = "0"))]
  data.table::setorder(matr_count, -N)
  for (colnum in seq(along = matr_names)) {
    colname = matr_names[colnum]
    renaming2 = matr_count[duplicated(get(colname)) == F]
    colname_v2 = paste0(colname, "_v2")
    renaming2[, `:=`((colname_v2), paste(col1_factor_level, 
                                         get(colname)))]
    matr_count[, `:=`((colname_v2), renaming2[match_hk(matr_count[, 
                                                                  get(colname)], renaming2[, get(colname)]), get(colname_v2)])]
  }
  alluvial2 = function(..., freq, col = "gray", border = 0, 
                       layer, hide = FALSE, alpha = 0.5, gap.width = 0.05, 
                       xw = 0.1, cw = 0.1, blocks = TRUE, ordering = NULL, 
                       axis_labels = NULL, cex = par("cex"), cex.axis = par("cex.axis"), 
                       leadingLetters2remove = 6 + 1) {
    p <- data.frame(..., freq = freq, col, alpha, border, 
                    hide, stringsAsFactors = FALSE)
    np <- ncol(p) - 5
    if (!is.null(ordering)) {
      stopifnot(is.list(ordering))
      if (length(ordering) != np) 
        stop("'ordering' argument should have ", np, 
             " components, has ", length(ordering))
    }
    n <- nrow(p)
    if (missing(layer)) {
      layer <- 1:n
    }
    p$layer <- layer
    d <- p[, 1:np, drop = FALSE]
    p <- p[, -c(1:np), drop = FALSE]
    p$freq <- with(p, freq/sum(freq))
    col <- col2rgb(p$col, alpha = TRUE)
    if (!identical(alpha, FALSE)) {
      col["alpha", ] <- p$alpha * 256
    }
    p$col <- apply(col, 2, function(x) do.call(rgb, c(as.list(x), 
                                                      maxColorValue = 256)))
    isch <- sapply(d, is.character)
    d[isch] <- lapply(d[isch], as.factor)
    if (length(blocks) == 1) {
      blocks <- if (!is.na(as.logical(blocks))) {
        rep(blocks, np)
      }
      else if (blocks == "bookends") {
        c(TRUE, rep(FALSE, np - 2), TRUE)
      }
    }
    if (is.null(axis_labels)) {
      axis_labels <- names(d)
    }
    else {
      if (length(axis_labels) != ncol(d)) 
        stop("`axis_labels` should have length ", names(d), 
             ", has ", length(axis_labels))
    }
    getp <- function(i, d, f, w = gap.width) {
      a <- c(i, (1:ncol(d))[-i])
      if (is.null(ordering[[i]])) {
        o <- do.call(order, d[a])
      }
      else {
        d2 <- d
        d2[1] <- ordering[[i]]
        o <- do.call(order, d2[a])
      }
      x <- c(0, cumsum(f[o])) * (1 - w)
      x <- cbind(x[-length(x)], x[-1])
      gap <- cumsum(c(0L, diff(as.numeric(d[o, i])) != 
                        0))
      mx <- max(gap)
      if (mx == 0) 
        mx <- 1
      gap <- gap/mx * w
      (x + gap)[order(o), ]
    }
    dd <- lapply(seq_along(d), getp, d = d, f = p$freq)
    rval <- list(endpoints = dd)
    op <- par(mar = mymargins)
    plot(NULL, type = "n", xlim = c(1 - cw, np + cw), ylim = c(0, 
                                                               1), xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", 
         xlab = "", ylab = "", frame = FALSE)
    ind <- which(!p$hide)[rev(order(p[!p$hide, ]$layer))]
    for (i in ind) {
      for (j in 1:(np - 1)) {
        xspline(c(j, j, j + xw, j + 1 - xw, j + 1, j + 
                    1, j + 1 - xw, j + xw, j) + rep(c(cw, -cw, 
                                                      cw), c(3, 4, 2)), c(dd[[j]][i, c(1, 2, 2)], 
                                                                          rev(dd[[j + 1]][i, c(1, 1, 2, 2)]), dd[[j]][i, 
                                                                                                                      c(1, 1)]), shape = c(0, 0, 1, 1, 0, 0, 1, 
                                                                                                                                           1, 0, 0), open = FALSE, col = p$col[i], border = p$border[i])
      }
    }
    for (j in seq_along(dd)) {
      ax <- lapply(split(dd[[j]], d[, j]), range)
      names(ax) = stringr::str_sub(string = names(ax), 
                                   start = leadingLetters2remove, end = stringr::str_length(names(ax)))
      if (blocks[j]) {
        for (k in seq_along(ax)) {
          rect(j - cw, ax[[k]][1], j + cw, ax[[k]][2])
        }
      }
      else {
        for (i in ind) {
          x <- j + c(-1, 1) * cw
          y <- t(dd[[j]][c(i, i), ])
          w <- xw * (x[2] - x[1])
          xspline(x = c(x[1], x[1], x[1] + w, x[2] - 
                          w, x[2], x[2], x[2] - w, x[1] + w, x[1]), 
                  y = c(y[c(1, 2, 2), 1], y[c(2, 2, 1, 1), 
                                            2], y[c(1, 1), 1]), shape = c(0, 0, 1, 
                                                                          1, 0, 0, 1, 1, 0, 0), open = FALSE, col = p$col[i], 
                  border = p$border[i])
        }
      }
      for (k in seq_along(ax)) {
        text(j, mean(ax[[k]]), labels = names(ax)[k], 
             cex = cex)
      }
    }
    axis(1, at = rep(c(-cw, cw), ncol(d)) + rep(seq_along(d), 
                                                each = 2), line = 0.5, col = "white", col.ticks = "black", 
         labels = FALSE)
    axis(1, at = seq_along(d), tick = FALSE, labels = axis_labels, 
         cex.axis = cex.axis)
    par(op)
    invisible(rval)
  }
  data.table::setorder(matr_count, -N)
  for (i in matr_names) matr_count[, `:=`((i), NULL)]
  data.table::setnames(matr_count, paste0(matr_names, "_v2"), 
                       matr_names)
  if (colorlines == F) 
    alluvial2(matr_count[, matr_names, with = F], freq = matr_count$N, 
              col = rainbow(matr_count[, uniqueN(get(spalte4color))])[factor(matr_count[, 
                                                                                        get(spalte4color)]) %>% as.numeric], gap.width = gap.width, 
              border = "grey90", alpha = alpha, blocks = showblocks, 
              cw = x_axis_size, ...)
  if (colorlines == T) 
    alluvial2(matr_count[, matr_names, with = F], freq = matr_count$N, 
              col = rainbow(matr_count[, uniqueN(get(spalte4color))])[factor(matr_count[, 
                                                                                        get(spalte4color)]) %>% as.numeric], gap.width = gap.width, 
              border = scales::alpha(rainbow(matr_count[, uniqueN(get(spalte4color))]), 
                                     alpha = min(0.25, alpha - 0.55))[factor(matr_count[, 
                                                                                        get(spalte4color)]) %>% as.numeric], alpha = alpha, 
              blocks = showblocks, cw = x_axis_size, ...)
}


# # LOAD seurat objekt laden----
seurat_fn = here("results/s603_2_seurat_lunge_merged_21-12-11.RDS")
stopifnot(file.exists(seurat_fn))
seurat_fn



Errint_harmony =  readRDS(seurat_fn)
Errint_harmony@active.ident = factor(Errint_harmony$SCT_snn_res.0.3)

Errint_harmony$celltype <- paste("C", Errint_harmony$SCT_snn_res.0.3)
Errint_harmony@active.ident = factor(Errint_harmony$celltype)

Errint_harmony <- RenameIdents (Errint_harmony, "C 0" = "T Cells", 'C 1' = "Endothelial", 'C 2' = "Alveolar Macrophages", 'C 3' = "Macrophages/ Monocytes", 'C 4' = "NK Cells", 'C 5' = "AT2", 'C 6' = "Endothelial", 'C 7' = "Fibroblasts", 'C 8' = "B Cells", 'C 9' = "Perivascular", 'C 10' = "AT1", 'C 11' = "Macrophages/ Monocytes", 'C 12' = "Mast Cells", 'C 13' = "Ciliated Cells", 'C 14' = "Club Cells", 'C 15' = "Dendritic Cells", 'C 16' = "Macrophages/ Monocytes", 'C 17' = "ly Endothelial", 'C 18' = "Proliferating NK/T", "C 19" = "Proliferating AM", "C 20" = "Neutrophils", "C 21" = "Endothelial", "C 22" = "Fibroblasts")



reihenfolge = c('Alveolar Macrophages', 'Proliferating AM', 'Macrophages/ Monocytes', "Neutrophils", 'Dendritic Cells', 'Mast Cells', 'T Cells', 'NK Cells', 'Proliferating NK/T', 'B Cells', 'AT1', 'AT2', 'Ciliated Cells', 'Club Cells', 'Endothelial', 'ly Endothelial', 'Fibroblasts', 'Perivascular')

NogpaletteReihe <-  c("#CB769E", "#DE639A", "#A85C85", "#0081AF", "#4F6D7A", "#7C6A0A", "#368F8B", "#246A73", "#5CC1BC", "#62C370", "#F7C548", "#F97E44", "#FB3640", "#B7245C", "#0D3B66", "#3E2F5B", "#B2675E", "#644536")

venn2(Errint_harmony@active.ident, reihenfolge)
Errint_harmony@active.ident <- factor(x = Errint_harmony@active.ident, levels = reihenfolge)




p_harmony = UMAPPlot(Errint_harmony, label = T, cols = c(NogpaletteReihe), label.size = 3.0) + NoLegend()  + ggtitle("Integration via Harmony")+ xlab("UMAP 1") + ylab("UMAP 2") + scale_color_manual(drop=FALSE, values = NogpaletteReihe)
p_harmony


# alternative laden
Errint_seurat = readRDS(here("results/s630_1_seurat_lunge_merged_noMITOadjSEURATAnchors_22-01-02.RDS"))

Errint_seurat_anno = Errint_seurat@meta.data %>% as.data.table()

Errint_seurat_anno[,uniqueN(integrated_snn_res.0.3)]
Errint_seurat_anno[,uniqueN(integrated_snn_res.0.4)]

Errint_seurat_anno[, .(prozNeutro = sum(rev1_celltype =="Neutrophils")), integrated_snn_res.0.4][order(prozNeutro)]
Errint_seurat_anno[integrated_snn_res.0.4 == 4,.N,rev1_celltype ]


venn2(Errint_seurat$rev1_celltype_completed, reihenfolge)
Errint_seurat@active.ident = factor(Errint_seurat$rev1_celltype_completed, levels  = reihenfolge)

levels(Errint_harmony@active.ident)
levels(Errint_seurat@active.ident)

p_seurat = UMAPPlot(Errint_seurat, label = T, cols = c(NogpaletteReihe), label.size = 3.0) + NoLegend() + ggtitle("Integration via Seurat integrateData()") + xlab("UMAP 1") + ylab("UMAP 2") + scale_color_manual(drop=FALSE, values = NogpaletteReihe)
p_seurat

require(patchwork)
p_harmony + p_seurat

Errint_seurat_anno_plot = Errint_seurat_anno[,.N, .(harmony = rev1_celltype, integrateData = rev1_celltype_completed)]


require(ggsankey)
require(ggrepel)

require(dplyr)
Errint_seurat_anno_plot2 = ggsankey::make_long(Errint_seurat_anno[,.(harmony = rev1_celltype, integrateData = rev1_celltype_completed)], harmony, integrateData)

Errint_seurat_anno_plot2$next_node = factor(Errint_seurat_anno_plot2$next_node, levels = reihenfolge)


p_sankey = ggplot(Errint_seurat_anno_plot2, aes(x = x, next_x = next_x, 
                     node = node, next_node = next_node, 
                     fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.6, node.color = alpha("grey55", 0.2), alpha = 0.4) +
  geom_sankey_text(size = 4, color = "black") +
  scale_fill_manual(values = NogpaletteReihe) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
                    plot.title = element_text(hjust = .5),
        axis.text.x = element_text(color = "black"))
p_sankey


require(patchwork)
 ((p_harmony / p_seurat) | p_sankey) + patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")") &  theme(plot.tag = element_text(face = 'bold', size = 14))

pdf(here("results/s631_2_supplFig_vgl_UMAP_harmony_seuratAnchor.pdf"), 9.7, 11.5)
((p_harmony / p_seurat) | p_sankey) + patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")") &  theme(plot.tag = element_text(face = 'bold', size = 14))
dev.off()



finalizeSkript()
