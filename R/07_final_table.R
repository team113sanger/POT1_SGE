## =============================================================================
## Stage 6 - the final annotated variant tables
##
## Takes the stage-5 variant table and attaches the four annotation sources
## that cannot be computed from the raw reads, then applies the final column
## selection, naming and ordering.
##
##   SpliceAI scores and GRCh38-normalised coordinates
##   EVE and CADD scores
##   FoldX/MutFunc ddG predictions
##   plasmid libraries for targetons 3_1a and 14b, which were sequenced
##     separately from the rest
##
## The first three are shipped in data/external/; see data/README.md for how
## each was produced.
##
## Output
##   supplementary_data_2  one row per annotated oligo
##   supplementary_data_1  one row per variant
## =============================================================================

## One variant cannot be lifted to GRCh38: it is an indel that also carries a
## PAM edit, so its coordinates are ambiguous. It has no SpliceAI annotation and
## is therefore absent from the final tables. The join below asserts that it is
## the only row lost there, so a future change to the annotation file that
## silently dropped other variants would stop the pipeline rather than pass.
UNMAPPABLE_ID <-
  "ENST00000357628.8.ENSG00000128513.16_chr7:124892332_124892347_clinvar_20201107_rc"

## Replicates excluded at the DESeq2 stage are still present in the annotation
## frame. Blank them so the reported counts show only what was actually used.
blank_excluded_replicates <- function(df) {
  r3 <- stringr::str_detect(df$source_targetons, "9_1")
  r1 <- stringr::str_detect(df$source_targetons, "14")
  d7 <- stringr::str_detect(df$source_targetons, "10")
  for (tp in params$timepoints) {
    df[[paste0(tp, "R3")]][r3] <- NA
    df[[paste0(tp, "R1")]][r1] <- NA
  }
  df$D7R3[d7] <- NA
  df
}

## Final column order and names. Note that "functional_classificaiton" is
## misspelt; the spelling is retained deliberately so that column names stay
## stable for anything already reading these tables.
FINAL_RENAME <- c(
  pos = "fixed_pos", ref = "fixed_ref", alt = "fixed_alt",
  functional_classificaiton = "Cluster", exon = "EXON", intron = "INTRON",
  vep_consequence_slim = "slim_consequence", vep_consequence = "Consequence",
  Ref_AA = "ref_aa", Alt_AA = "alt_aa",
  target_region = "REGION", target_region_data_source = "source_targetons",
  variation_id = "VariationID", clinical_significance = "CLIN_SIG.y",
  clinical_significance_slim = "ClinVar_CLNSIG_slim",
  LFC_D4_D7 = "processed_LFC_D4_D7", LFC_D4_D10 = "processed_LFC_D4_D10",
  LFC_D4_D14 = "processed_LFC_D4_D14", LFC_D4_D21 = "processed_LFC_D4_D21",
  LFC_continuous = "processed_LFC_continuous",
  BH_FDR_D4_D7 = "processed_BH_FDR_D4_D7", BH_FDR_D4_D10 = "processed_BH_FDR_D4_D10",
  BH_FDR_D4_D14 = "processed_BH_FDR_D4_D14", BH_FDR_D4_D21 = "processed_BH_FDR_D4_D21",
  BH_FDR_continuous = "processed_BH_FDR_continuous",
  SE_D4_D7 = "SE_bind_D4_D7", SE_D4_D10 = "SE_bind_D4_D10",
  SE_D4_D14 = "SE_bind_D4_D14", SE_D4_D21 = "SE_bind_D4_D21",
  SE_continuous = "SE_bind_continuous",
  z_D4_D7 = "processed_Z_D4_D7", z_D4_D10 = "processed_Z_D4_D10",
  z_D4_D14 = "processed_Z_D4_D14", z_D4_D21 = "processed_Z_D4_D21",
  z_continuous = "processed_Z_continuous")

FINAL_COLUMNS <- c(
  "HGVSc","HGVSp","chrom_pos_ref_alt","pos","ref","alt",
  "functional_classificaiton","functional_score","exon","intron",
  "vep_consequence_slim","cDNA_position","CDS_position","Protein_position",
  "Amino_acids","Ref_AA","Alt_AA","Codons","gene_id","transcript_id","mseq","id",
  "target_region","target_region_data_source","functional_score_below_threshold",
  "primer_flag","pam_codon","mutator","mave_nt","vep_consequence","DOMAINS",
  "HGVS_OFFSET",
  paste0("LFC_",    LFC_TAGS), paste0("BH_FDR_", LFC_TAGS),
  paste0("SE_",     LFC_TAGS), paste0("z_",      LFC_TAGS),
  paste0("two_tailed_p_", LFC_TAGS),
  "C1_z","C2_z","uncertainty", sample_columns(), "PLASMID",
  "SIFT","PolyPhen","EVE_SCORE","EVE_CLASS","CADD_PHRED","CADD_RAW","ddG",
  "DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL",
  "REF_DS_AG","REF_DS_AL","REF_DS_DG","REF_DS_DL",
  "ALT_DS_AG","ALT_DS_AL","ALT_DS_DG","ALT_DS_DL",
  "AF","gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF",
  "gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF","gnomAD_FLAG",
  "dbSNP","variation_id","is_in_clinvar","clinvar_onestar","clinical_significance",
  "clinical_significance_slim","is_in_gnomad","present_accession",
  "present_accession_filter")

## How the tables are laid out: mutator (snv first, then the factor order used
## used throughout), then descending genomic position, then oligo id so
## nothing depends on the order rows happen to arrive in.
canonical_order <- function(df) {
  order(factor(df$mutator, levels = MUTATOR_LEVELS), -df$pos, df$id, method = "radix")
}

## How a representative row is chosen when a variant appears more than once.
## The same keys, with one extra preference inserted: a measurement combined
## across both tiled targetons beats a single-targeton one, since it pools two
## independent observations. Ordering on mutator first is what makes a variant
## reachable as a plain SNV be reported as one, rather than as a recoded or
## deleted codon.
selection_order <- function(df) {
  order(factor(df$mutator, levels = MUTATOR_LEVELS), -df$pos,
        !grepl("|", df$target_region_data_source, fixed = TRUE),
        df$id, method = "radix")
}

MUTATOR_LEVELS <- c("snv","1del","custom","2del0","2del1","inframe","stop","ala","snvre")

## Columns added after the selection step, so excluded from it.
DERIVED_COLUMNS <- c("chrom_pos_ref_alt", "functional_score_below_threshold", "ddG")

## The upstream column name behind each final column name.
source_name <- function(final) {
  ifelse(final %in% names(FINAL_RENAME), FINAL_RENAME[final], final)
}

build_final_table <- function(variants) {

  ## ---- external: SpliceAI and GRCh38-normalised coordinates --------------
  sp <- utils::read.csv(file.path(paths$external, "spliceai_fixed_coords.csv.gz"),
                        check.names = FALSE)
  key_v <- paste(variants$mseq, variants$id)
  key_s <- paste(sp$mseq, sp$id)
  dropped <- setdiff(key_v, key_s)
  if (length(dropped) != 1L || variants$id[match(dropped, key_v)] != UNMAPPABLE_ID)
    stop("expected exactly the unmappable variant to be absent from the SpliceAI ",
         "annotation, got ", length(dropped), " unmatched row(s)")

  df <- merge(variants, sp,
              by = c("mseq", "id"), all.x = FALSE)

  ## ---- blank the excluded replicates --------------------------------------
  df <- blank_excluded_replicates(df)

  ## ---- external: EVE and CADD ---------------------------------------------
  eve <- utils::read.csv(file.path(paths$external, "eve_cadd_vep.csv.gz"))
  eve <- eve[!is.na(eve$EVE_SCORE) & eve$EVE_SCORE != "-", ]
  eve$EVE_SCORE <- as.numeric(eve$EVE_SCORE)
  if (anyDuplicated(eve$HGVSc)) stop("EVE annotation has duplicate HGVSc keys")
  df <- merge(df, eve, by = "HGVSc", all.x = TRUE)

  ## ---- reduce to the reported columns, then rename -------------------------
  ## Selecting first matters: the VaLiAnT metafile carries its own ref and new
  ## allele columns, which would collide with fixed_ref/fixed_alt -> ref/alt.
  keep <- source_name(setdiff(FINAL_COLUMNS, DERIVED_COLUMNS))
  missing <- setdiff(keep, names(df))
  if (length(missing)) stop("source columns absent: ", paste(missing, collapse = ", "))
  df <- df[, keep]
  for (new in names(FINAL_RENAME)) {
    old <- FINAL_RENAME[[new]]
    if (old %in% names(df)) names(df)[names(df) == old] <- new
  }

  ## ---- derived columns -------------------------------------------------------
  df$chrom_pos_ref_alt <- paste0("7_", df$pos, "_", df$ref, "_", df$alt)

  ## VaLiAnT targeton numbering -> exon numbering
  df$target_region <- unname(REGION_TO_EXON[as.character(df$target_region)])
  df$target_region_data_source <-
    unname(SOURCE_TARGETON_TO_EXON[as.character(df$target_region_data_source)])
  if (anyNA(df$target_region) || anyNA(df$target_region_data_source))
    stop("unmapped targeton label when renumbering to exon numbers")

  ## Depleted variants below the functional-score threshold.
  ## Left NA for unchanged and enriched variants, where it has no meaning.
  df$functional_score_below_threshold <-
    ifelse(df$functional_classificaiton %in% c("unchanged", "enriched"), NA,
           ifelse(df$functional_score < params$score_threshold, "T", "F"))

  ## Asn->Ala is "N/A" in single-letter code, which R and Excel both read as
  ## missing. Written back explicitly.
  df$Amino_acids[df$Codons %in% c("AAT/GCC", "AAC/GCC")] <- "N/A"

  ## ---- plasmid counts for the two separately sequenced libraries -----------
  p3 <- utils::read.delim(file.path(paths$plasmid, "pot1_exon_3_1_a_lib.counts.tsv.gz"))
  p3 <- data.frame(id = p3$NAME, PLASMID_3 = p3$reads_pot1_exon_3_1_a_lib)
  p14 <- utils::read.delim(file.path(paths$plasmid, "POT1_14_B_PLASMID_LIB.lib_counts.tsv.gz"))
  p14 <- data.frame(id = p14$NAME, PLASMID_14 = p14$COUNT)
  if (anyDuplicated(p3$id) || anyDuplicated(p14$id))
    stop("plasmid library count files contain duplicate oligo names")

  df <- merge(df, p3,  by = "id", all.x = TRUE)
  df <- merge(df, p14, by = "id", all.x = TRUE)
  df$PLASMID <- dplyr::case_when(
    !is.na(df$PLASMID)                                           ~ as.numeric(df$PLASMID),
    stringr::str_detect(df$target_region_data_source, "8_1|8-1") ~ as.numeric(df$PLASMID_3),
    stringr::str_detect(df$target_region_data_source, "19")      ~ as.numeric(df$PLASMID_14),
    TRUE                                                         ~ NA_real_)
  df$PLASMID_3 <- NULL; df$PLASMID_14 <- NULL

  ## ---- external: FoldX/MutFunc ddG -----------------------------------------
  dg <- utils::read.csv(file.path(paths$external, "ddG_mutfunc.csv.gz"))
  dg$mapping <- paste0(dg$aa_wt, "_", dg$uniprot_pos, "_", dg$aa_mt)
  if (anyDuplicated(dg$mapping))
    stop("ddG table has duplicate residue mappings; the join would duplicate rows")
  df$mapping <- paste0(df$Ref_AA, "_", df$Protein_position, "_", df$Alt_AA)
  df <- merge(df, dg[, c("mapping", "ddG")], by = "mapping", all.x = TRUE)
  df$mapping <- NULL

  ## ---- final row order -------------------------------------------------------
  df$mutator <- factor(df$mutator, levels = MUTATOR_LEVELS)
  df <- df[canonical_order(df), , drop = FALSE]
  rownames(df) <- NULL

  missing <- setdiff(FINAL_COLUMNS, names(df))
  if (length(missing)) stop("final table is missing: ", paste(missing, collapse = ", "))
  df[, FINAL_COLUMNS]
}

## Supplementary Data 1: one row per variant.
##
## Reduced to one row per mutant sequence, then to one row per HGVSc, keeping
## the first row in the canonical order above. Eleven coding changes are
## reachable twice: eight are pairs of adjacent single-base deletions that
## produce the same sequence (the higher position is kept), and three are the
## PAM-edited and unedited forms of the same SNV in the 3_1/3_2 overlap (the
## measurement combined across both targetons is kept).
unique_final_table <- function(final) {
  u <- final[selection_order(final), , drop = FALSE]
  u <- u[!duplicated(u$mseq), , drop = FALSE]
  u <- u[!duplicated(u$HGVSc), , drop = FALSE]
  u <- u[canonical_order(u), , drop = FALSE]
  rownames(u) <- NULL
  u
}

## Write a workbook: title in row 1, column header in row 2, and a second sheet
## documenting every column. The titles and the notes sheet are editorial text
## held in data/metadata/ rather than generated here.
write_supplementary <- function(data, n, file) {
  ttl   <- readLines(file.path(paths$metadata, sprintf("supp_data_%d_titles.txt", n)))
  notes <- utils::read.csv(file.path(paths$metadata, sprintf("supp_data_%d_notes.csv", n)),
                           check.names = FALSE)
  wb <- openxlsx::createWorkbook()
  s1 <- sprintf("Supplementary Data %d", n)
  openxlsx::addWorksheet(wb, s1)
  openxlsx::writeData(wb, s1, ttl[1], startRow = 1, colNames = FALSE)
  openxlsx::writeData(wb, s1, data,   startRow = 2)
  openxlsx::addWorksheet(wb, ttl[3])
  openxlsx::writeData(wb, ttl[3], ttl[2], startRow = 1, colNames = FALSE)
  openxlsx::writeData(wb, ttl[3], CLINICAL_USE_NOTICE, startRow = 2, colNames = FALSE)
  openxlsx::writeData(wb, ttl[3], notes, startRow = 4)
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}
