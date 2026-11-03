## =============================================================================
## Optional QC plots, one PDF set per library. Off by default: they add several
## minutes per library and none of them feed the results.
## =============================================================================

scree_plot <- function(rld, ntop = 500) {
  m  <- SummarizedExperiment::assay(rld)
  rv <- matrixStats::rowVars(m)
  sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- stats::prcomp(t(m[sel, ]))
  pv  <- pca$sdev^2 / sum(pca$sdev^2)
  xx  <- graphics::barplot(round(pv, 2) * 100, names.arg = seq_along(pv),
                           xlab = "PC", ylab = "% Variance", ylim = c(0, 100),
                           main = "Scree Plot")
  graphics::text(x = xx, y = round(pv, 4) * 100, label = round(pv, 4) * 100,
                 pos = 3, cex = 0.8)
}

ggpairs_ext <- function(data, mapping, pts = list(), smt = list(), ...) {
  ggplot2::ggplot(data = data, mapping = mapping, ...) +
    do.call(ggplot2::geom_point, pts) + do.call(ggplot2::geom_smooth, smt)
}

write_qc_plots <- function(dds, rld, lib) {
  id <- lib$lib_id
  f  <- function(x) file.path(paths$qc, paste0(id, x))

  grDevices::pdf(f("_scree.pdf"));      scree_plot(rld);                   grDevices::dev.off()
  grDevices::pdf(f("_dispersion.pdf")); DESeq2::plotDispEsts(dds, ylim = c(1e-4, 2e1)); grDevices::dev.off()

  d <- as.matrix(stats::dist(t(SummarizedExperiment::assay(rld))))
  grDevices::pdf(f("_heatmap.pdf"))
  pheatmap::pheatmap(d, trace = "none",
    col = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255))
  grDevices::dev.off()

  grDevices::pdf(f("_scatter_matrix.pdf"), width = 8, height = 8)
  for (tp in params$timepoints) {
    sub <- as.data.frame(SummarizedExperiment::assay(rld))
    sub <- sub[, grepl(paste0("^", tp, "R"), colnames(sub)), drop = FALSE]
    if (ncol(sub) < 2) next
    print(GGally::ggpairs(sub,
      lower = list(continuous = GGally::wrap(ggpairs_ext,
        pts = list(size = 0.4, colour = "black"),
        smt = list(method = "lm", se = FALSE, linewidth = 0.2, colour = "blue"))),
      title = paste0(tp, " - regularised-log read count")))
  }
  grDevices::dev.off()

  pca <- DESeq2::plotPCA(rld, intgroup = c("condition", "type"), returnData = TRUE)
  pv  <- round(100 * attr(pca, "percentVar"))
  p <- ggplot2::ggplot(pca, ggplot2::aes(PC1, PC2, color = condition, shape = type)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::xlab(paste0("PC1: ", pv[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", pv[2], "% variance")) +
    ggplot2::coord_fixed() + ggplot2::theme_classic()
  ggplot2::ggsave(f("_PCA.pdf"), p, height = 6, width = 8)

  utils::write.table(as.data.frame(SummarizedExperiment::assay(rld)), sep = "\t",
                     file = f("_rld.txt"), col.names = NA)
  utils::write.table(DESeq2::counts(dds, normalized = TRUE), sep = "\t",
                     file = f("_norm_count.txt"), col.names = NA)
}
