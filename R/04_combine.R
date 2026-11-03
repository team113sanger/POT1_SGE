## =============================================================================
## Stage 4 - combine libraries into one variant-level table
##
## Only the 15 libraries marked in_final contribute. Where an exon was covered
## by both guide libraries only one was carried forward, so the only variants
## measured twice are those in the overlap of the tiled targeton pairs
## (3_1/3_2 and 9_1/9_2). Those two measurements are combined by inverse-
## variance weighting; everything else passes through unchanged.
##
## VaLiAnT targeton label -> REGION. Tiled pairs share a REGION, which is what
## makes them group together for the weighted combination below.
## =============================================================================

EXON_GROUP_TO_REGION <- c("1"="2","2"="3","3_1"="4","3_2"="4","4"="5","5"="6",
                          "6"="7","8"="9","9_1"="10","9_2"="10","10"="11",
                          "11"="12","12"="13","13"="14","14"="15")

combine_libraries <- function(deseq_results, manifest = library_manifest) {

  final_libs <- manifest$lib_id[manifest$in_final]
  tiled      <- intersect(final_libs, TILED_LIBS)

  parts <- lapply(final_libs, function(l) {
    d <- deseq_results[[l]]$out
    d$process <- if (l %in% tiled) "tiled" else "non_tiled"
    d
  })
  total <- plyr::rbind.fill(parts)

  ## Which guide library's PAM edit, if any, touches this mutant sequence.
  ## A wins over B where both are present.
  pam <- total %>%
    dplyr::group_by(mseq) %>%
    dplyr::summarize(
      dup_mseq = dplyr::n(),
      PAM_status_full_arrayed = paste(pam_mut_sgrna_id, collapse = ","),
      .groups = "drop") %>%
    dplyr::mutate(PAM_status_simplified = dplyr::case_when(
      stringr::str_detect(PAM_status_full_arrayed, "_a") ~ "A",
      stringr::str_detect(PAM_status_full_arrayed, "_b") ~ "B",
      TRUE ~ "-")) %>%
    dplyr::select(mseq, dup_mseq, PAM_status_simplified)

  tbl <- merge(total, pam, by = "mseq", all.x = TRUE) %>%
    dplyr::arrange(EXON_GROUP, id, SG) %>%
    dplyr::rename(Variant_Sources = SG,
                  PAM_status = pam_mut_sgrna_id,
                  PAM_status_simplified_collated = PAM_status_simplified,
                  dup_mseq_collated = dup_mseq) %>%
    dplyr::rename_with(~ sub("^uncombined_two_tailed_p_", "two_tailed_p_", .x),
                       dplyr::starts_with("uncombined_two_tailed_p_")) %>%
    dplyr::mutate(
      Variant_duplication = NA,
      dataset_process = "not_combined",
      mseq_a = ifelse(Variant_Sources == "A", mseq, "NA"),
      mseq_b = ifelse(Variant_Sources == "B", mseq, "NA"),
      REGION = unname(EXON_GROUP_TO_REGION[EXON_GROUP]))

  if (anyNA(tbl$REGION))
    stop("EXON_GROUP values with no REGION mapping: ",
         paste(unique(tbl$EXON_GROUP[is.na(tbl$REGION)]), collapse = ", "))

  tbl$mseq <- NULL

  ## ---- collapse to one row per (variant id, region) ----------------------
  collapsed <- tbl %>%
    dplyr::group_by(id, REGION) %>%
    dplyr::summarize(
      n_targetons = dplyr::n(),
      source_dataset   = paste(dataset_process, collapse = "|"),
      source_targetons = paste(EXON_GROUP, collapse = "|"),
      source_libraries = paste(Variant_Sources, collapse = "|"),
      variant_n_instances_between_libraries = paste(Variant_duplication, collapse = "|"),
      PAM_status_simplified_collated = paste(PAM_status_simplified_collated, collapse = "|"),
      PAM_status = paste(PAM_status, collapse = "|"),
      dup_mseq_collated = paste(dup_mseq_collated, collapse = "|"),
      arrayed_mseq_a = paste(mseq_a, collapse = "|"),
      arrayed_mseq_b = paste(mseq_b, collapse = "|"),
      dplyr::across(dplyr::all_of(paste0("lfcSE_",   LFC_TAGS)), ~ paste(.x, collapse = ",")),
      dplyr::across(dplyr::all_of(paste0("adj_lfc_", LFC_TAGS)), ~ paste(.x, collapse = ",")),
      .groups = "drop")

  ## split the collapsed "primary,secondary" strings back into numeric pairs
  for (tg in LFC_TAGS) {
    for (stem in c("lfcSE", "adj_lfc")) {
      col <- paste0(stem, "_", tg)
      sp  <- do.call(rbind, strsplit(paste0(collapsed[[col]], ","), ",", fixed = TRUE))
      collapsed[[paste0("primary_tg_",   stem, "_", tg)]] <- as.numeric(sp[, 1])
      collapsed[[paste0("secondary_tg_", stem, "_", tg)]] <- suppressWarnings(as.numeric(sp[, 2]))
      collapsed[[col]] <- NULL
    }
  }

  ## ---- inverse-variance weighted combination ------------------------------
  single <- collapsed$n_targetons == 1L
  for (tg in LFC_TAGS) {
    se1 <- collapsed[[paste0("primary_tg_lfcSE_",   tg)]]
    se2 <- collapsed[[paste0("secondary_tg_lfcSE_", tg)]]
    l1  <- collapsed[[paste0("primary_tg_adj_lfc_",   tg)]]
    l2  <- collapsed[[paste0("secondary_tg_adj_lfc_", tg)]]

    w1 <- ifelse(single, NA_real_, 1 / se1^2)
    w2 <- ifelse(single, NA_real_, 1 / se2^2)
    wsum <- w1 + w2                       # NA if either is missing, by design

    collapsed[[paste0("SE_bind_", tg)]] <-
      ifelse(single, se1, wsum^(-0.5))
    collapsed[[paste0("processed_LFC_", tg)]] <-
      ifelse(single, l1, (w1 * l1 + w2 * l2) / wsum)
  }

  ## ---- pick the annotation source and a representative mseq ---------------
  ## Library A is preferred wherever the variant was seen in A.
  collapsed <- collapsed %>%
    dplyr::mutate(required_annotation_source = dplyr::case_when(
      stringr::str_detect(source_libraries, "A") ~ "A",
      stringr::str_detect(source_libraries, "B") ~ "B",
      TRUE ~ "-"))

  split_first <- function(x) vapply(strsplit(x, "|", fixed = TRUE), `[`, character(1), 1L)
  collapsed$mseq_a_tg1 <- split_first(collapsed$arrayed_mseq_a)
  collapsed$mseq_b_tg1 <- split_first(collapsed$arrayed_mseq_b)
  collapsed$mseq_for_filter <- ifelse(collapsed$required_annotation_source == "B",
                                      collapsed$mseq_b_tg1, collapsed$mseq_a_tg1)

  if (any(collapsed$mseq_for_filter == "NA"))
    stop("some variants resolved to a placeholder mseq - annotation source logic is wrong")

  ## one row per mutant sequence, so that the FDR correction below counts
  ## each independent observation once
  collapsed <- collapsed[!duplicated(collapsed$mseq_for_filter), , drop = FALSE]

  ## ---- final statistics ---------------------------------------------------
  for (tg in LFC_TAGS) {
    z <- collapsed[[paste0("processed_LFC_", tg)]] / collapsed[[paste0("SE_bind_", tg)]]
    p <- stats::pnorm(abs(z), lower.tail = FALSE) * 2
    collapsed[[paste0("processed_Z_", tg)]]      <- z
    collapsed[[paste0("two_tailed_p_", tg)]]     <- p
    collapsed[[paste0("processed_BH_FDR_", tg)]] <- stats::p.adjust(p, method = "BH")
  }

  collapsed$functional_classification <- dplyr::case_when(
    collapsed$processed_BH_FDR_continuous >= params$fdr_threshold ~ "unchanged",
    collapsed$processed_LFC_continuous < 0 ~ "depleted",
    collapsed$processed_LFC_continuous > 0 ~ "enriched",
    TRUE ~ NA_character_)

  if (anyNA(collapsed$functional_classification))
    warning(sum(is.na(collapsed$functional_classification)),
            " variants could not be classified")

  collapsed
}
