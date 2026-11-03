## =============================================================================
## POT1 SGE analysis — configuration
##
## Everything the pipeline needs to know about where files live, which
## libraries were sequenced, and what the analysis constants are. Nothing is
## configured anywhere else.
##
## All paths are relative to the repository root, so the pipeline runs wherever
## the repository is cloned. Run it from the repository root, or let
## run_pipeline.R resolve the root for you.
## =============================================================================

## Resolved by run_pipeline.R; defaults to the working directory so the file can
## also be sourced interactively from the repository root.
if (!exists("REPO_ROOT")) REPO_ROOT <- normalizePath(".", mustWork = TRUE)

paths <- list(
  ## --- experimental data -------------------------------------------------
  counts         = file.path(REPO_ROOT, "data", "counts"),
  plasmid        = file.path(REPO_ROOT, "data", "plasmid"),

  ## --- library design and variant annotation -----------------------------
  library_design = file.path(REPO_ROOT, "data", "library_design"),
  vep            = file.path(REPO_ROOT, "data", "vep"),

  ## --- annotation from outside this experiment ---------------------------
  external       = file.path(REPO_ROOT, "data", "external"),
  metadata       = file.path(REPO_ROOT, "data", "metadata"),

  ## --- outputs -------------------------------------------------------------
  results        = file.path(REPO_ROOT, "results"),
  per_library    = file.path(REPO_ROOT, "results", "per_library"),
  qc             = file.path(REPO_ROOT, "results", "qc"),
  logs           = file.path(REPO_ROOT, "results", "logs")
)

for (d in c(paths$results, paths$per_library, paths$qc, paths$logs))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)


## =============================================================================
## Library manifest
##
## Eighteen libraries were sequenced: nine against guide library A and nine
## against guide library B, each over five timepoints in triplicate.
##
## region_start / region_end
##   GRCh38 coordinates of the targeton. They are not decoration: the VaLiAnT
##   design and VEP annotation filenames are derived from them further down, so
##   the manifest is the single place a library is described.
##
## plasmid_mode
##   How the plasmid library was sequenced, which also determines how the SGE
##   counts are laid out on disk:
##     "merged"    one pre-merged plasmid file; 15 count files, one per sample
##     "two_lane"  plasmid sequenced over two lanes, and the SGE counts split
##                 across the same two lanes: 30 count files, summed in pairs
##     "none"      no plasmid library was sequenced for this targeton, so
##                 PLASMID is recorded as NA
##
## in_final
##   Whether the library contributes to the combined variant table. Where a
##   targeton was covered by both guide libraries only one was carried forward
##   (targeton 3_2 -> A, 4 -> A, 6 -> B). The other three are still processed
##   through DESeq2, so their QC output is available, but they are excluded
##   downstream to avoid counting the same variant twice from the same region.
##
## drop_reps
##   Sample columns excluded for that library after quality control. Three
##   libraries lost a replicate: the experimental designs in build_design() are
##   generated from this column, so an exclusion is declared once and flows
##   through the count matrix, the normalisation matrix and the model.
## =============================================================================

library_manifest <- tibble::tribble(
  ~exon,  ~sg, ~region_start, ~region_end, ~plasmid_mode, ~in_final, ~drop_reps,
  "2",    "a", 124870880L, 124871127L, "merged",   TRUE,  character(0),
  "3_1",  "a", 124863272L, 124863520L, "none",     TRUE,  character(0),
  "3_2",  "a", 124863468L, 124863717L, "two_lane", TRUE,  character(0),
  "4",    "a", 124858890L, 124859138L, "two_lane", TRUE,  character(0),
  "6",    "a", 124851781L, 124852030L, "two_lane", FALSE, character(0),
  "8",    "a", 124842757L, 124843006L, "merged",   TRUE,  character(0),
  "10",   "a", 124835214L, 124835463L, "merged",   TRUE,  "D7R3",
  "11",   "a", 124829216L, 124829465L, "merged",   TRUE,  character(0),
  "12",   "a", 124827112L, 124827359L, "merged",   TRUE,  character(0),
  "1",    "b", 124892164L, 124892411L, "merged",   TRUE,  character(0),
  "3_2",  "b", 124863468L, 124863717L, "two_lane", FALSE, character(0),
  "4",    "b", 124858890L, 124859138L, "two_lane", FALSE, character(0),
  "5",    "b", 124852932L, 124853177L, "merged",   TRUE,  character(0),
  "6",    "b", 124851781L, 124852030L, "two_lane", TRUE,  character(0),
  "9_1",  "b", 124840935L, 124841182L, "merged",   TRUE,  c("D4R3","D7R3","D10R3","D14R3","D21R3"),
  "9_2",  "b", 124840984L, 124841232L, "merged",   TRUE,  character(0),
  "13",   "b", 124825226L, 124825474L, "merged",   TRUE,  character(0),
  "14",   "b", 124823919L, 124824168L, "none",     TRUE,  c("D4R1","D7R1","D10R1","D14R1","D21R1")
)

## Targetons 3_1/3_2 and 9_1/9_2 tile the same exon and overlap, so a variant
## in the overlap is measured twice. Those pairs of measurements are combined by
## inverse-variance weighting in stage 4. Every other library stands alone.
TILED_LIBS <- c("E3_1_SGA", "E3_2_SGA", "E9_1_SGB", "E9_2_SGB")

library_manifest <- dplyr::mutate(
  library_manifest,
  lib_id    = paste0("E", exon, "_SG", toupper(sg)),
  ## VaLiAnT writes targeton 3_1 as "w3.1"; the count directories use "3_1"
  valiant_w = paste0("w", gsub("_", ".", exon), "_", sg),
  region    = paste0("chr7_", region_start, "_", region_end, "_minus_sgRNA_", valiant_w),
  meta_file = file.path(paths$library_design, paste0(region, "_meta.csv.gz")),
  vep_file  = file.path(paths$vep, paste0("VEP_OUTPUT_", region, "_pam.tsv.gz")),
  count_dir = file.path(paths$counts, lib_id)
)


## =============================================================================
## Analysis constants
## =============================================================================

params <- list(
  ## POT1 MANE Select transcript. VEP reports every transcript; all annotation
  ## is restricted to this one.
  transcript        = "ENST00000357628.8",

  timepoints        = c("D4", "D7", "D10", "D14", "D21"),
  replicates        = c("R1", "R2", "R3"),

  ## Days elapsed since D4. Used as the covariate in the continuous model,
  ## whose slope is the reported functional score.
  timepoint_days    = c(D4 = 0, D7 = 3, D10 = 6, D14 = 10, D21 = 17),

  ## Oligos with fewer than this many reads summed across the 15 SGE samples
  ## are dropped before modelling.
  min_row_count     = 10,

  ## DESeq2 ignores oligos at or below this total when estimating size factors
  ## from the control set.
  size_factor_minrowsum = 10,

  ## Variants expected to be functionally neutral. Size factors are estimated
  ## from these alone, and log fold changes are centred on their median, so a
  ## genuine genome-wide depletion is not normalised away.
  neutral_consequences = c("synonymous_variant", "intron_variant"),

  ## Depleted / enriched / unchanged, on the continuous model.
  fdr_threshold     = 0.01,

  ## Gaussian mixture splitting the depleted variants into strong and weak.
  gmm_model         = "EEV",     # ellipsoidal, equal volume, varying orientation
  gmm_k             = 2L,
  gmm_seed          = 0L,
  gmm_vars          = c("processed_LFC_D4_D21", "processed_LFC_D4_D7"),

  ## Functional-score cut-off applied within the depleted class.
  score_threshold   = -0.03359,

  ## 5' primer-binding sites. Variants inside these windows are flagged, since
  ## a variant under the primer can distort that oligo's read count.
  primer_windows    = list(c(124892384, 124892388), c(124870905, 124870909),
                           c(124859113, 124859117), c(124825250, 124825254)),

  ## File formats: pycroquet writes two comment lines before the header, and
  ## VEP --tab writes 92 lines of ## metadata.
  pycroquet_skip    = 2L,
  vep_skip          = 92L
)

## Carried into the Notes sheet of each output workbook, so that it travels with
## the file if the workbook is downloaded and circulated on its own.
CLINICAL_USE_NOTICE <- paste(
  "These functional scores are research data from a study that has not yet completed",
  "peer review. They must not be used for clinical diagnosis, variant interpretation,",
  "patient management, or any other clinical purpose.")

## The four discrete contrasts plus the continuous slope. Used as a suffix
## throughout the statistics columns.
LFC_TAGS <- c("D4_D7", "D4_D10", "D4_D14", "D4_D21", "continuous")

## Three mutant sequences in the 3_1/3_2 overlap are generated by a delins
## mutator and duplicate a sequence already present from a different oligo.
## They survive deduplication because the oligo ids differ, so they are removed
## by name. All three are classified unchanged.
MANUAL_ID_EXCLUSIONS <- paste0(
  "ENST00000357628.8.ENSG00000128513.16_chr7:124863497_124863499_TAC>",
  c("GAC", "CAC", "AAC"), "_snvre_rc")

## VaLiAnT numbers targetons from the start of the designed region; the
## reported tables number exons. These map one onto the other.
REGION_TO_EXON <- c("2"="6","3"="7","4"="8","5"="9","6"="10","7"="11",
                    "9"="13","10"="14","11"="15","12"="16","13"="17",
                    "14"="18","15"="19")

SOURCE_TARGETON_TO_EXON <- c(
  "1"="6", "2"="7", "3_1"="8_1", "3_1|3_2"="8-1|8_2", "3_2"="8_2",
  "4"="9", "5"="10", "6"="11", "8"="13", "9_1"="14_1", "9_1|9_2"="14_1|14_2",
  "9_2"="14_2", "10"="15", "11"="16", "12"="17", "13"="18", "14"="19")
