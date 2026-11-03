## =============================================================================
## Stage 3 - DESeq2 per library
##
## Two models are fitted per library against the same size factors:
##   discrete   ~ condition  with D4 as reference, giving D4->D7/D10/D14/D21
##   continuous ~ condition  with condition = days since D4, giving one slope
##                           ("the SGE rate")
##
## Size factors are estimated from the neutral variants only (synonymous and
## intronic), then forced onto the full model, so that a genuine global
## depletion is not normalised away.
##
## Log fold changes are then median-centred on those same neutral variants and
## the p-values recomputed, because the median-centring shifts the estimate
## without changing its standard error.
##
## Dispersion trend: DESeq2 fits a parametric trend where it can and falls back
## to a local regression where it cannot. That choice is made by an iterative
## Gamma GLM started at a fixed point, c(0.1, 1), which errors out if any
## coefficient goes non-positive on the way. For a library whose asymptotic
## dispersion is near zero the default start overshoots, the fit throws, and
## DESeq2 silently substitutes a different model - a decision that can flip on a
## ~1e-9 change in the inputs, i.e. on compiled code and hardware rather than on
## the data. fit_dispersion_trend() below removes that sensitivity.
## =============================================================================

## Experimental designs are derived from the manifest rather than read from
## per-library CSVs. These reproduce annot_a.csv, annot_a_9_1.csv,
## annot_a_14.csv and annot_a_10.csv exactly.
build_design <- function(lib) {
  keep <- setdiff(sample_columns(), lib$drop_reps[[1]])
  tp   <- sub("R[123]$", "", keep)
  rep  <- sub("^D[0-9]+",  "", keep)
  list(
    samples   = keep,
    discrete  = data.frame(condition = tp,
                           type = rep, sample = "sg", row.names = NULL),
    continuous = data.frame(condition = as.numeric(params$timepoint_days[tp]),
                            type = rep, sample = "sg", row.names = NULL)
  )
}

## Count matrix keyed on mseq. Oligos sharing an mseq are collapsed to the
## first occurrence: different VaLiAnT mutators can generate the same mutant
## sequence, but the physical library contains only one copy of it.
as_count_matrix <- function(df, samples) {
  df <- df[!duplicated(df$mseq), , drop = FALSE]
  m  <- as.matrix(df[, samples, drop = FALSE])
  rownames(m) <- df$mseq
  m
}

## The neutral set used for normalisation.
neutral_subset <- function(df) {
  df <- df[!duplicated(df$mseq), , drop = FALSE]
  df[df$Consequence %in% params$neutral_consequences, , drop = FALSE]
}

## Size factors from control oligos only.
## DESeq2's parametric dispersion fit, with the starting value exposed.
## Identical to DESeq2:::parametricDispersionFit otherwise.
parametric_fit_from <- function(means, disps, start) {
  coefs <- start
  for (iter in 0:10) {
    residuals <- disps / (coefs[1] + coefs[2] / means)
    good <- which(residuals > 1e-4 & residuals < 15)
    suppressWarnings(
      fit <- stats::glm(disps[good] ~ I(1 / means[good]),
                        family = stats::Gamma(link = "identity"), start = coefs))
    old <- coefs
    coefs <- stats::coefficients(fit)
    if (!all(coefs > 0)) return(NULL)
    if (sum(log(coefs / old)^2) < 1e-6 && fit$converged) {
      names(coefs) <- c("asymptDisp", "extraPois")
      return(coefs)
    }
  }
  NULL
}

## Try the parametric trend from DESeq2's own starting value first, then from
## progressively smaller asymptotic dispersions, and only fall back to a local
## regression if none of them converge. The solution is the same wherever it is
## reached - only the path to it is fragile - so this removes the hardware
## sensitivity without changing the model being fitted.
##
## Eight of the eighteen libraries fit parametrically under this rule and ten
## fall back to a local regression. The rule is deterministic: the fit type
## follows from the data, not from the machine the pipeline runs on.
PARAMETRIC_STARTS <- list(c(0.1, 1), c(0.01, 1), c(0.001, 2), c(1e-4, 2.5), c(1e-5, 3))

fit_dispersion_trend <- function(dds, quiet = TRUE) {
  dds <- DESeq2::estimateDispersionsGeneEst(dds, quiet = quiet)
  md  <- S4Vectors::mcols(dds)
  use <- which(!md$allZero & !is.na(md$dispGeneEst) & md$dispGeneEst > 100 * 1e-8)

  coefs <- NULL
  for (st in PARAMETRIC_STARTS) {
    coefs <- parametric_fit_from(md$baseMean[use], md$dispGeneEst[use], st)
    if (!is.null(coefs)) break
  }

  if (is.null(coefs)) {
    dds <- DESeq2::estimateDispersionsFit(dds, fitType = "local", quiet = quiet)
  } else {
    dispFn <- function(q) coefs[1] + coefs[2] / q
    attr(dispFn, "coefficients") <- coefs
    attr(dispFn, "fitType") <- "parametric"
    DESeq2::dispersionFunction(dds) <- dispFn
  }
  DESeq2::estimateDispersionsMAP(dds, quiet = quiet)
}

## DESeq() with the dispersion trend fitted as above.
deseq_fit <- function(dds, quiet = TRUE) {
  dds <- fit_dispersion_trend(dds, quiet = quiet)
  DESeq2::nbinomWaldTest(dds, quiet = quiet)
}

estimate_control_size_factors <- function(countData, colData, design,
                                          minRowSum = params$size_factor_minrowsum,
                                          ref = NULL) {
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = colData,
                                        design = stats::as.formula(design))
  dds <- dds[rowSums(DESeq2::counts(dds)) > minRowSum, ]
  if (!is.null(ref)) dds$condition <- stats::relevel(dds$condition, ref = ref)
  DESeq2::sizeFactors(DESeq2::estimateSizeFactors(dds))
}

append_suffix <- function(df, suffix) {
  colnames(df) <- paste0(colnames(df), suffix)
  df
}

run_deseq_library <- function(annotated, lib, qc_plots = FALSE) {
  des <- build_design(lib)
  samples <- des$samples

  count_mat   <- as_count_matrix(annotated, samples)
  neutral_mat <- as_count_matrix(neutral_subset(annotated), samples)

  ## ---- discrete model ----------------------------------------------------
  sf <- estimate_control_size_factors(neutral_mat, des$discrete,
                                      design = "~ condition", ref = "D4")

  dds <- DESeq2::DESeqDataSetFromMatrix(count_mat, des$discrete, design = ~ condition)
  dds$condition <- factor(dds$condition, levels = params$timepoints)
  dds$condition <- stats::relevel(dds$condition, ref = "D4")
  dds <- DESeq2::estimateSizeFactors(dds)
  DESeq2::sizeFactors(dds) <- sf
  dds <- deseq_fit(dds)
  rld <- DESeq2::rlog(dds)

  ## per-variant z-scores of the regularised-log counts
  z <- as.data.frame(t(scale(t(as.matrix(SummarizedExperiment::assay(rld))))))
  colnames(z) <- paste0(colnames(z), "_z_score")
  z$mseq <- rownames(z); rownames(z) <- NULL
  z <- z[, c("mseq", setdiff(names(z), "mseq"))]

  ## No LFC shrinkage: the "raw" estimates are used throughout. apeglm would
  ## be the RNA-seq default but shrinks real SGE effects toward zero.
  contrasts <- paste0("condition_", c("D7", "D10", "D14", "D21"), "_vs_D4")
  tw <- DEGreport::degComps(dds, combs = "condition", contrast = as.list(contrasts),
                            alpha = 0.05, skip = FALSE, type = "normal",
                            pairs = FALSE, fdr = "default")

  suffixes <- c("_D4_D7", "_D4_D10", "_D4_D14", "_D4_D21")
  parts <- Map(function(i, sfx) {
    d <- as.data.frame(DEGreport::deg(tw[[i]], "raw"))
    d <- append_suffix(d, sfx)
    tibble::rownames_to_column(d, var = "mseq")
  }, seq_along(suffixes), suffixes)

  summary_tbl <- Reduce(function(a, b) dplyr::left_join(a, b, by = "mseq"), parts)
  summary_tbl <- summary_tbl[, !names(summary_tbl) %in%
                               c("baseMean_D4_D10", "baseMean_D4_D14", "baseMean_D4_D21")]
  names(summary_tbl)[names(summary_tbl) == "baseMean_D4_D7"] <- "baseMean"
  summary_tbl <- merge(summary_tbl, z, by = "mseq", all.x = TRUE)

  ## ---- continuous model --------------------------------------------------
  sf_c <- estimate_control_size_factors(neutral_mat, des$continuous,
                                        design = "~ condition")
  ddc <- DESeq2::DESeqDataSetFromMatrix(count_mat, des$continuous, design = ~ condition)
  ddc <- DESeq2::estimateSizeFactors(ddc)
  DESeq2::sizeFactors(ddc) <- sf_c
  ddc <- deseq_fit(ddc)
  twc <- DEGreport::degComps(ddc, combs = "condition", alpha = 0.05, skip = FALSE,
                             type = "normal", pairs = FALSE, fdr = "default")
  rate <- tibble::rownames_to_column(as.data.frame(twc[[1]]), var = "mseq")
  rate <- append_suffix(rate, "_continuous")
  names(rate)[names(rate) == "mseq_continuous"] <- "mseq"
  summary_tbl <- merge(summary_tbl, rate, by = "mseq", all.x = TRUE)

  ## ---- median-centre on the neutral set ----------------------------------
  lfc_cols <- c(paste0("log2FoldChange", suffixes), "log2FoldChange_continuous")
  se_cols  <- c(paste0("lfcSE",          suffixes), "lfcSE_continuous")
  tags     <- c(sub("^_", "", suffixes), "continuous")

  neutral_mseq <- rownames(neutral_mat)
  neutral_rows <- summary_tbl[match(neutral_mseq, summary_tbl$mseq), , drop = FALSE]
  medians <- vapply(lfc_cols, function(cc) stats::median(neutral_rows[[cc]]), numeric(1))

  for (k in seq_along(tags)) {
    tg <- tags[k]
    summary_tbl[[paste0("median_", tg)]]    <- medians[[lfc_cols[k]]]
    summary_tbl[[paste0("adj_lfc_", tg)]]   <- summary_tbl[[lfc_cols[k]]] - medians[[lfc_cols[k]]]
    summary_tbl[[paste0("adj_score_", tg)]] <- summary_tbl[[paste0("adj_lfc_", tg)]] /
                                                summary_tbl[[se_cols[k]]]
    p <- stats::pnorm(abs(summary_tbl[[paste0("adj_score_", tg)]]), lower.tail = FALSE) * 2
    summary_tbl[[paste0("uncombined_two_tailed_p_", tg)]] <- p
    summary_tbl[[paste0("uncombined_BH_FDR_", tg)]] <- stats::p.adjust(p, method = "BH")
  }

  ## ---- assemble outputs ---------------------------------------------------
  out <- summary_tbl
  out$SG         <- toupper(lib$sg)
  out$IDENT      <- lib$lib_id
  out$EXON_GROUP <- lib$exon
  out <- merge(out, annotated[, c("id", "mseq", "pam_mut_sgrna_id")],
               by = "mseq", all.x = TRUE)
  out$combined_name <- paste0(out$id, "_", out$EXON_GROUP)

  out_full <- merge(out, annotated, by = "id", all.x = TRUE)
  names(out_full)[names(out_full) == "pam_mut_sgrna_id.x"] <- "pam_mut_sgrna_id"
  names(out_full)[names(out_full) == "mseq.x"] <- "mseq"
  out_full$pam_mut_sgrna_id.y <- NULL
  out_full$mseq.y <- NULL

  if (qc_plots) write_qc_plots(dds, rld, lib)

  list(out = out, out_full = out_full, dds = dds, rld = rld,
       size_factors = sf, medians = medians)
}

run_deseq_all <- function(annotated, manifest = library_manifest, qc_plots = FALSE) {
  out <- lapply(seq_len(nrow(manifest)), function(i) {
    lib <- as.list(manifest[i, ])
    message("  DESeq2: ", lib$lib_id,
            if (length(lib$drop_reps[[1]])) paste0("  (excluding ",
              paste(lib$drop_reps[[1]], collapse = ","), ")") else "")
    r <- run_deseq_library(annotated[[lib$lib_id]], lib, qc_plots = qc_plots)
    utils::write.csv(r$out, file.path(paths$per_library, paste0(lib$lib_id, "_OUT.csv")),
                     row.names = FALSE)
    utils::write.csv(r$out_full, file.path(paths$per_library,
                     paste0(lib$lib_id, "_OUT_full_anot.csv")), row.names = FALSE)
    r
  })
  stats::setNames(out, manifest$lib_id)
}
