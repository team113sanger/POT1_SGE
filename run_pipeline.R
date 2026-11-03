#!/usr/bin/env Rscript
## =============================================================================
## POT1 saturation genome editing — analysis pipeline
##
##   Rscript run_pipeline.R [--qc-plots] [--stage N]
##
## Reads the raw per-library read counts in data/ and produces the variant-level
## functional score tables in results/. Nothing outside this repository is read
## or written, and no path needs editing: everything resolves relative to the
## repository root.
##
##   --stage N    stop after stage N (default 6, the whole pipeline)
##   --qc-plots   additionally write the per-library QC plots to results/qc/
##                (several minutes per library; none of it feeds the results)
##
## Stages
##   1  read counts                       results/per_library/*_count_frame.csv
##   2  attach library design and VEP annotation, filter
##   3  model each library with DESeq2     results/per_library/*_OUT.csv
##   4  combine libraries                  results/pot1_combined_stats.csv
##   5  expand annotation, ClinVar, depletion classes
##                                         results/pot1_variant_table.csv
##   6  external annotation, final layout  results/supplementary_data_{1,2}.xlsx
##
## Runtime is roughly four minutes end to end on a laptop.
## =============================================================================

t_start <- Sys.time()

args      <- commandArgs(trailingOnly = TRUE)
qc_plots  <- "--qc-plots" %in% args
stage_arg <- if ("--stage" %in% args) as.integer(args[which(args == "--stage") + 1L]) else 6L

## Resolve the repository root from this script's own location, so the pipeline
## can be launched from anywhere.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
REPO_ROOT <- if (length(this_file)) dirname(normalizePath(this_file[1])) else
             normalizePath(".", mustWork = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(stringr); library(tidyr)
  library(DESeq2); library(DEGreport); library(mclust); library(openxlsx)
  library(ggplot2); library(GGally); library(pheatmap); library(RColorBrewer)
})

source(file.path(REPO_ROOT, "config.R"))
for (f in sprintf("0%d_%s.R", 1:7, c("counts", "annotate", "deseq", "combine",
                                     "classify", "qc_plots", "final_table")))
  source(file.path(REPO_ROOT, "R", f))

message("POT1 SGE analysis pipeline")
message("  repository : ", REPO_ROOT)
message("  libraries  : ", nrow(library_manifest),
        " (", sum(library_manifest$in_final), " contribute to the variant table)")
message("  QC plots   : ", qc_plots)

## ---- 1. counts ---------------------------------------------------------------
message("\n[1/6] reading counts")
count_frames <- build_all_count_frames()

## ---- 2. library design, annotation and filtering -----------------------------
message("\n[2/6] attaching library design and VEP annotation, filtering")
annotated <- annotate_all(count_frames)
fsum <- filter_summary(annotated)
write.csv(fsum, file.path(paths$results, "filter_summary.csv"), row.names = FALSE)
message("  oligos designed: ", sum(fsum$oligos_in),
        " | removed below ", params$min_row_count, " reads: ", sum(fsum$removed_low_count),
        " | removed in constant regions: ", sum(fsum$removed_const_region),
        " | retained: ", sum(fsum$after_const_region))
if (stage_arg < 3) quit(save = "no")

## ---- 3. per-library models ----------------------------------------------------
message("\n[3/6] fitting DESeq2 models")
deseq_results <- run_deseq_all(annotated, qc_plots = qc_plots)
if (stage_arg < 4) quit(save = "no")

## ---- 4. combine across libraries ----------------------------------------------
message("\n[4/6] combining libraries")
combined <- combine_libraries(deseq_results)
write.csv(combined, file.path(paths$results, "pot1_combined_stats.csv"), row.names = FALSE)
message("  unique variants: ", nrow(combined))
print(table(combined$functional_classification))
if (stage_arg < 5) quit(save = "no")

## ---- 5. full annotation and depletion classes ---------------------------------
message("\n[5/6] expanding annotation, ClinVar, depletion classes")
expanded <- expand_and_annotate(combined, annotated)
variants <- classify_depletion(expanded)
write.csv(variants, file.path(paths$results, "pot1_variant_table.csv"), row.names = FALSE)
uniq <- variants[!duplicated(variants$mseq), ]
message("  annotated rows: ", nrow(variants), " | unique mutant sequences: ", nrow(uniq))
print(table(uniq$Cluster))
if (stage_arg < 6) quit(save = "no")

## ---- 6. external annotation and final tables ----------------------------------
message("\n[6/6] external annotation and final tables")
final  <- build_final_table(variants)
final1 <- unique_final_table(final)
message("  Supplementary Data 2: ", nrow(final),  " rows x ", ncol(final),  " columns")
message("  Supplementary Data 1: ", nrow(final1), " rows x ", ncol(final1), " columns")

write.csv(final,  file.path(paths$results, "supplementary_data_2.csv"), row.names = FALSE)
write.csv(final1, file.path(paths$results, "supplementary_data_1.csv"), row.names = FALSE)
write_supplementary(final,  2, file.path(paths$results, "supplementary_data_2.xlsx"))
write_supplementary(final1, 1, file.path(paths$results, "supplementary_data_1.xlsx"))

## Every intermediate object, for anyone who wants to pick the analysis apart
## without re-running it.
saveRDS(list(count_frames = count_frames, annotated = annotated,
             deseq_results = lapply(deseq_results, `[`, c("out", "out_full",
                                                          "size_factors", "medians")),
             combined = combined, variants = variants, final = final),
        file.path(paths$results, "pipeline_objects.rds"))

## Exact package versions for this run.
writeLines(capture.output(sessionInfo()), file.path(paths$logs, "sessionInfo.txt"))

message("\ndone in ", round(difftime(Sys.time(), t_start, units = "mins"), 1), " min")
message("results written to ", paths$results)
