#!/usr/bin/env Rscript
## Install everything run_pipeline.R needs.
##
##   Rscript install_dependencies.R
##
## DESeq2 and DEGreport come from Bioconductor, the rest from CRAN.

cran <- c("dplyr", "tibble", "stringr", "tidyr", "purrr", "plyr", "mclust",
          "openxlsx", "ggplot2", "GGally", "pheatmap", "RColorBrewer",
          "cowplot", "ggpubr", "ggridges", "matrixStats", "locfit", "BiocManager")
bioc <- c("DESeq2", "DEGreport", "SummarizedExperiment", "S4Vectors")

missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE)

cat("\nInstalled versions:\n")
for (p in c(cran, bioc))
  cat(sprintf("  %-22s %s\n", p,
      tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")))
