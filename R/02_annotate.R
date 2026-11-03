## =============================================================================
## Stage 2 - attach the VaLiAnT design metadata and the VEP annotation,
##           then apply the count filters
## =============================================================================

## VaLiAnT metafile: one row per designed oligo, keyed on oligo_name.
read_meta <- function(lib) {
  if (!file.exists(lib$meta_file)) stop("missing metafile: ", lib$meta_file)
  m <- utils::read.csv(lib$meta_file, header = TRUE, comment.char = "#")
  names(m)[names(m) == "oligo_name"] <- "id"
  m
}

## VEP --tab output, restricted to the POT1 MANE transcript.
read_vep <- function(lib) {
  if (!file.exists(lib$vep_file)) stop("missing VEP output: ", lib$vep_file)
  v <- utils::read.table(lib$vep_file, skip = params$vep_skip, header = TRUE,
                         comment.char = "")
  names(v)[names(v) == "X.Uploaded_variation"] <- "id"
  v[v$Feature == params$transcript, , drop = FALSE]
}

## Join counts to design metadata. An inner join: oligos present in the count
## file but not in the design are dropped.
attach_meta <- function(counts, lib) {
  merge(counts, read_meta(lib), by = "id", all.x = FALSE)
}

## Two filters, applied in order:
##   1. total reads across the 15 SGE samples >= 10
##   2. the variant does not fall in a constant (non-targeted) region
## PLASMID is excluded from the row sum - it is a library QC measure, not a
## timepoint, and two libraries have no plasmid data at all.
filter_counts <- function(df, min_row = params$min_row_count) {
  sc <- sample_columns()
  df$rowsum <- rowSums(df[, sc])
  n0 <- nrow(df)
  df <- df[df$rowsum >= min_row, , drop = FALSE]
  n1 <- nrow(df)
  df <- df[df$vcf_var_in_const == 0 | is.na(df$vcf_var_in_const), , drop = FALSE]
  n2 <- nrow(df)
  attr(df, "filter_counts") <- c(input = n0, after_min_count = n1,
                                 after_const_region = n2)
  df
}

## Counts + design + VEP. all.x = TRUE: an oligo with no VEP record on the MANE
## transcript is kept, with NA annotation.
annotate_library <- function(counts, lib) {
  master   <- attach_meta(counts, lib)
  filtered <- filter_counts(master)
  fstat    <- attr(filtered, "filter_counts")
  annotated <- merge(filtered, read_vep(lib), by = "id", all.x = TRUE)
  attr(annotated, "filter_counts") <- fstat
  annotated
}

annotate_all <- function(count_frames, manifest = library_manifest) {
  out <- lapply(seq_len(nrow(manifest)), function(i) {
    lib <- as.list(manifest[i, ])
    message("  annotate: ", lib$lib_id)
    annotate_library(count_frames[[lib$lib_id]], lib)
  })
  stats::setNames(out, manifest$lib_id)
}

## Per-library filtering summary, written to results/filter_summary.csv so the
## number of oligos lost at each step is a reported output rather than
## something that has to be recovered from the console.
filter_summary <- function(annotated) {
  do.call(rbind, lapply(names(annotated), function(n) {
    a <- attr(annotated[[n]], "filter_counts")
    data.frame(lib_id = n, oligos_in = a[["input"]],
               after_min_count = a[["after_min_count"]],
               after_const_region = a[["after_const_region"]],
               removed_low_count = a[["input"]] - a[["after_min_count"]],
               removed_const_region = a[["after_min_count"]] - a[["after_const_region"]],
               row.names = NULL)
  }))
}
