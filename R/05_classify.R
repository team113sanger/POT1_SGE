## =============================================================================
## Stage 5 - expand back to full annotation, attach ClinVar, and split the
##           depleted variants into strong and weak by Gaussian mixture
## =============================================================================

## Stage 4 reduced each variant to one representative mutant sequence. Here it
## is expanded again so that every annotated oligo id carries the statistics.
build_meta_vep <- function(annotated, manifest = library_manifest) {
  final <- manifest[manifest$in_final, ]
  parts <- lapply(seq_len(nrow(final)), function(i) {
    lib <- as.list(final[i, ])
    d <- annotated[[lib$lib_id]]
    d$SG <- toupper(lib$sg)
    d$exon_group <- lib$exon
    d$required_annotation_id <- paste0(d$id, "_", d$SG)
    d
  })
  plyr::rbind.fill(parts)
}

## ClinVar. Two downloads taken on the same date: the full set of POT1 records,
## and the subset whose review status is one star or better. Both are shipped
## in data/external/ so the annotation is fixed rather than moving under the
## analysis; see data/README.md for the release date and query.
read_clinvar <- function() {
  strip_paren <- function(x) gsub("\\s*\\([^\\)]+\\)", "", as.character(x))
  to_enst <- function(x) sub("NM_015450.3", params$transcript, x, fixed = TRUE)

  rd <- function(f) utils::read.delim(f, header = TRUE, sep = "\t", quote = "",
                                      comment.char = "", check.names = TRUE)

  all <- rd(file.path(paths$external, "clinvar_pot1_all_20240312.txt.gz"))
  all$CLIN_SIG <- strip_paren(all$Germline.classification)
  all$HGVSc    <- to_enst(strip_paren(all$Name))
  all <- all[grepl(params$transcript, all$HGVSc, fixed = TRUE),
             c("VariationID", "CLIN_SIG", "HGVSc")]

  one <- rd(file.path(paths$external, "clinvar_pot1_onestar_20240312.txt.gz"))
  one$clinvar_onestar <- "Y"
  one$HGVSc <- to_enst(strip_paren(one$Name))
  one <- one[grepl(params$transcript, one$HGVSc, fixed = TRUE),
             c("clinvar_onestar", "HGVSc")]

  m <- merge(all, one, by = "HGVSc", all.x = TRUE)
  m$clinvar_onestar <- ifelse(is.na(m$clinvar_onestar), "N", "Y")
  m
}

## VEP reports a long Consequence string; collapse it to one label per variant.
## Order matters: later rules override earlier ones.
slim_consequence <- function(df) {
  cons <- df$Consequence
  mut  <- df$mutator
  s <- rep(NA_character_, nrow(df))
  hit <- function(pat) stringr::str_detect(cons, pat) & !is.na(cons)

  s[hit("stop_gained")]             <- "stop_gained"
  s[hit("UTR_variant")]             <- "UTR"
  s[hit("frameshift_variant")]      <- "frameshift"
  s[stringr::str_detect(mut, "inframe") & !is.na(mut)] <- "codon_deletion"
  s[hit("intron_variant")]          <- "intron"
  s[hit("missense_variant")]        <- "missense"
  s[hit("synonymous_variant")]      <- "synonymous"
  s[hit("splice_acceptor_variant")] <- "splice_acceptor"
  s[hit("splice_donor_variant")]    <- "splice_donor"

  simple <- mut %in% c("snvre", "snv", "custom", "ala")
  s[hit("start_lost") & simple] <- "start_lost"
  s[hit("stop_lost")  & simple] <- "stop_lost"
  s[hit("stop_retained_variant") & mut %in% c("snvre", "snv", "custom")] <- "synonymous"
  s[hit("inframe_deletion")  & mut == "custom"] <- "clinical_inframe_deletion"
  s[hit("inframe_insertion") & mut == "custom"] <- "clinical_inframe_insertion"
  s
}

KEEP_COLUMNS <- c(
  "mseq","HGVSc","HGVSp","id.y","REGION","source_targetons","source_libraries",
  "PAM_status_simplified_collated","PAM_status",
  paste0("SE_bind_",       LFC_TAGS),
  paste0("processed_LFC_", LFC_TAGS),
  as.vector(rbind(paste0("processed_Z_",      LFC_TAGS),
                  paste0("two_tailed_p_",     LFC_TAGS),
                  paste0("processed_BH_FDR_", LFC_TAGS))),
  "functional_classification",
  sample_columns(), "PLASMID",
  "gene_id","transcript_id","ref_chr","ref_strand","ref_start","ref_end",
  "mut_position","ref","new","ref_aa","alt_aa","mut_type","mutator",
  "oligo_length","pam_mut_sgrna_id","mave_nt","vcf_var_in_const",
  "Consequence","cDNA_position","CDS_position","Protein_position","Amino_acids",
  "Codons","Existing_variation","SIFT","PolyPhen","EXON","INTRON","DOMAINS",
  "HGVS_OFFSET","AF","gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF",
  "gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF",
  "CLIN_SIG","gnomAD_FLAG","gnomAD_AF.1","ClinVar","ClinVar_CLNSIG",
  "ClinVar_CLNREVSTAT","dbSNP")

expand_and_annotate <- function(combined, annotated, manifest = library_manifest) {
  meta_vep <- build_meta_vep(annotated, manifest)

  stats_tbl <- dplyr::rename(combined, mseq = mseq_for_filter)
  expanded  <- merge(stats_tbl, meta_vep, by = "mseq")

  ## Duplicate mutant sequences in the 3_1/3_2 overlap that survive mseq-level
  ## deduplication because the PAM-edited and unedited forms differ. All are
  ## classified unchanged.
  expanded <- expanded[!expanded$id.y %in% MANUAL_ID_EXCLUSIONS, , drop = FALSE]

  missing <- setdiff(KEEP_COLUMNS, names(expanded))
  if (length(missing)) stop("expected columns absent: ", paste(missing, collapse = ", "))
  clean <- expanded[, KEEP_COLUMNS]

  clean$slim_consequence <- slim_consequence(clean)
  clean$functional_score <- clean$processed_LFC_continuous

  clinvar <- read_clinvar()
  clean <- merge(clean, clinvar, by = "HGVSc", all.x = TRUE)

  ## CLIN_SIG.x comes from the VEP annotation, CLIN_SIG.y from the ClinVar
  ## download; the download is the one used.
  cs <- clean$CLIN_SIG.y
  clean %>%
    dplyr::mutate(
      is_in_clinvar = ifelse(is.na(VariationID), "N", "Y"),
      is_in_gnomad  = ifelse(gnomAD_AF == "-", "N", "Y"),
      present_accession = dplyr::case_when(
        is_in_clinvar == "Y" & is_in_gnomad == "Y" ~ "ClinVar and gnomAD",
        is_in_clinvar == "Y"                       ~ "ClinVar only",
        is_in_gnomad  == "Y"                       ~ "gnomAD only",
        TRUE                                       ~ "Unobserved"),
      ClinVar_CLNSIG_slim = dplyr::case_when(
        is.na(cs)                                            ~ "Unobserved",
        cs %in% c("Pathogenic", "Likely pathogenic",
                  "Pathogenic/Likely pathogenic")            ~ "Pathogenic/Likely pathogenic",
        cs %in% c("Benign", "Likely benign",
                  "Benign/Likely benign")                    ~ "Benign/Likely benign",
        cs == "Uncertain significance"                       ~ "Uncertain significance",
        cs == "Conflicting classifications of pathogenicity" ~ "Conflicting interpretation",
        TRUE                                                 ~ NA_character_),
      present_accession_filter = paste0(present_accession, "_", mseq),
      present_accession = factor(present_accession,
        levels = c("ClinVar and gnomAD","ClinVar only","gnomAD only","Unobserved")),
      mutator = factor(mutator,
        levels = c("snv","1del","custom","2del0","2del1","inframe","stop","ala","snvre"))) %>%
    dplyr::arrange(present_accession) %>%
    dplyr::arrange(-mut_position) %>%
    dplyr::arrange(mutator) %>%
    dplyr::filter(!HGVSc %in% "-") %>%
    dplyr::rename(id = id.y)
}

## Depleted variants fall into two groups on the D4->D21 against D4->D7 plane:
## those that drop out early and those that drop out slowly. A two-component
## Gaussian mixture with ellipsoidal, equal-volume, varying-orientation
## covariance separates them. The seed is fixed in config.R so the fit is
## reproducible.
classify_depletion <- function(expanded) {
  uniq <- expanded[!duplicated(expanded$mseq), , drop = FALSE]

  depleted  <- uniq[uniq$functional_classification == "depleted", , drop = FALSE]
  enriched  <- uniq[uniq$functional_classification == "enriched", , drop = FALSE]
  unchanged <- uniq[uniq$functional_classification == "unchanged", , drop = FALSE]

  set.seed(params$gmm_seed)
  fit <- mclust::Mclust(depleted[, params$gmm_vars],
                        G = params$gmm_k, modelNames = params$gmm_model)
  message(sprintf("  GMM: %s, G=%d, BIC=%.3f, n=%d",
                  params$gmm_model, params$gmm_k, max(fit$BIC, na.rm = TRUE), nrow(depleted)))

  depleted <- cbind(depleted,
                    stats::setNames(as.data.frame(fit$z), c("C1_z", "C2_z")),
                    Cluster = fit$classification,
                    uncertainty = fit$uncertainty)
  enriched$Cluster  <- 3L
  unchanged$Cluster <- 4L

  rates <- plyr::rbind.fill(depleted, enriched, unchanged)
  rates$Cluster <- factor(rates$Cluster, levels = 1:4,
                          labels = c("strong depleted","weak depleted","enriched","unchanged"))
  rates$Cluster <- stats::relevel(rates$Cluster, "unchanged")

  out <- merge(expanded,
               rates[, c("mseq","C1_z","C2_z","Cluster","uncertainty")],
               by = "mseq")

  ## Flag variants within a 5 bp primer-binding window, where read counts may
  ## be distorted, and variants falling in a PAM-edited codon.
  in_window <- Reduce(`|`, lapply(params$primer_windows,
                                  function(w) out$mut_position >= w[1] & out$mut_position <= w[2]))
  out$primer_flag <- ifelse(in_window, "Y", "N")
  out$pam_codon   <- ifelse(stringr::str_detect(out$PAM_status, "sgRNA"), "Y", "N")

  out %>% dplyr::arrange(-mut_position) %>% dplyr::arrange(mutator)
}
