## =============================================================================
## Stage 1 - read pycroquet counts into one count matrix per library
##
## Each library contributes one count file per sample (or two, for libraries
## sequenced over two lanes). The timepoint and replicate are parsed out of the
## filename and the columns assembled by name, with an explicit check that all
## fifteen timepoint x replicate combinations are present. Naming the columns
## positionally from list.files() would also work today, but would fail
## silently if a file were ever renamed or missing.
## =============================================================================

sample_columns <- function() {
  as.vector(t(outer(params$timepoints, params$replicates, paste0)))
}

## Read one pycroquet query_to_library_counts file. Column 5 holds the count
## and is named after the sample; rename it to something predictable.
## Count files are stored gzipped and read directly; no decompression step is
## needed before running the pipeline.
read_pycroquet <- function(file, value_name = "count") {
  df <- utils::read.table(file, skip = params$pycroquet_skip, header = TRUE,
                          comment.char = "", fill = TRUE)
  df <- df[, c(1, 5)]
  names(df) <- c("id", value_name)
  df
}

## "2_a_D4_rep_1"  ->  D4R1
## "4a_D21_rep_3__44030_2_65" -> D21R3
parse_sample <- function(fname) {
  m <- regmatches(fname, regexec("_(D4|D7|D10|D14|D21)_rep_([123])", fname))[[1]]
  if (length(m) != 3L) return(NA_character_)
  paste0(m[2], "R", m[3])
}

## --- SGE counts ------------------------------------------------------------
read_sge_counts <- function(lib) {
  files <- list.files(lib$count_dir, pattern = "query_to_library_counts\\.tsv(\\.gz)?$",
                      full.names = TRUE)
  if (!length(files))
    stop("no count files in ", lib$count_dir)

  samples <- vapply(basename(files), parse_sample, character(1))
  if (anyNA(samples))
    stop(lib$lib_id, ": cannot parse timepoint/replicate from ",
         paste(basename(files)[is.na(samples)], collapse = ", "))

  frames <- Map(function(f, s) read_pycroquet(f, s), files, samples)

  ## Libraries sequenced over two lanes contribute two files per sample; sum
  ## them. Libraries sequenced on one lane contribute exactly one.
  expected_per_sample <- if (lib$plasmid_mode == "two_lane") 2L else 1L
  tab <- table(samples)
  if (!all(tab == expected_per_sample))
    stop(lib$lib_id, ": expected ", expected_per_sample, " file(s) per sample, got ",
         paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = " "))

  wanted <- sample_columns()
  if (!setequal(unique(samples), wanted))
    stop(lib$lib_id, ": sample set does not match the 15 expected timepoint x replicate combinations")

  counts <- Reduce(function(a, b) dplyr::full_join(a, b, by = "id"),
                   lapply(wanted, function(s) {
                     parts <- frames[samples == s]
                     out <- Reduce(function(a, b) dplyr::full_join(a, b, by = "id"),
                                   Map(function(p, i) stats::setNames(p, c("id", paste0(s, "_", i))),
                                       parts, seq_along(parts)))
                     out[[s]] <- rowSums(out[, -1, drop = FALSE], na.rm = TRUE)
                     out[, c("id", s)]
                   }))
  counts
}

## --- Plasmid counts --------------------------------------------------------
read_plasmid_counts <- function(lib) {
  if (lib$plasmid_mode == "none")
    return(NULL)

  if (lib$plasmid_mode == "merged") {
    pat <- sprintf("^Exon_%s_%s_merged\\.query_to_library_counts\\.tsv(\\.gz)?$",
                   lib$exon, lib$sg)
    f <- list.files(paths$plasmid, pattern = pat, full.names = TRUE)
    if (length(f) != 1L)
      stop(lib$lib_id, ": expected 1 merged plasmid file, found ", length(f))
    return(read_pycroquet(f, "PLASMID"))
  }

  ## two_lane: one file per lane, summed
  pat <- sprintf("^%s%s_lib__44030_[12]_[0-9]+\\.query_to_library_counts\\.tsv(\\.gz)?$",
                 lib$exon, lib$sg)
  f <- list.files(paths$plasmid, pattern = pat, full.names = TRUE)
  if (length(f) != 2L)
    stop(lib$lib_id, ": expected 2 plasmid lane files, found ", length(f))
  lanes <- Map(function(p, i) read_pycroquet(p, paste0("lane", i)), f, seq_along(f))
  pl <- Reduce(function(a, b) dplyr::full_join(a, b, by = "id"), lanes)
  pl$PLASMID <- rowSums(pl[, c("lane1", "lane2")], na.rm = TRUE)
  pl[, c("id", "PLASMID")]
}

## --- One count frame per library -------------------------------------------
build_count_frame <- function(lib) {
  counts  <- read_sge_counts(lib)
  plasmid <- read_plasmid_counts(lib)

  if (is.null(plasmid)) {
    counts$PLASMID <- NA_real_
  } else {
    counts <- dplyr::left_join(counts, plasmid, by = "id")
  }

  ## Oligos absent from a given file have no row rather than a zero.
  ## PLASMID is deliberately left as NA where no plasmid library exists.
  sc <- sample_columns()
  counts[sc][is.na(counts[sc])] <- 0
  if (lib$plasmid_mode != "none")
    counts$PLASMID[is.na(counts$PLASMID)] <- 0

  counts[, c("id", sc, "PLASMID")]
}

build_all_count_frames <- function(manifest = library_manifest, write = TRUE) {
  out <- lapply(seq_len(nrow(manifest)), function(i) {
    lib <- as.list(manifest[i, ])
    message("  counts: ", lib$lib_id)
    cf <- build_count_frame(lib)
    if (write)
      utils::write.csv(cf, file.path(paths$per_library,
                                     paste0(lib$lib_id, "_count_frame.csv")),
                       row.names = FALSE)
    cf
  })
  stats::setNames(out, manifest$lib_id)
}
