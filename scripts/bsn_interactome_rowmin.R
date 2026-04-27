source("R/load_fragpipe_utils.R")
source("R/sample_metadata_utils.R")
source("R/imputation_utils.R")
source("R/interactome_utils.R")
source("R/export_utils.R")

library(tidyverse)

df <- read_fragpipe_protein_table("data/abundance_protein_MD.tsv")

wt_pairs <- build_ip_pairs(
  replicate = c("TP652", "BN141", "BN368"),
  igg = c("TP652-WT-IgG", "BN141-WT-IgG", "BN368-WT-IgG"),
  ip = c("TP652-WT-BSN_IP", "BN141-WT-BSN_IP", "BN368-WT-BSN_IP"),
  condition = "WT"
)

homo_pairs <- build_ip_pairs(
  replicate = c("BN207", "BN206", "BN214"),
  igg = c("BN207-BSNKIHomo-IgG", "BN206-BSNKIHomo-IgG", "BN214-BSNKIHomo-IgG"),
  ip = c("BN207-BSNKIHomo-BSNIP", "BN206-BSNKIHomo-BSNIP", "BN214-BSNKIHomo-BSNIP"),
  condition = "Homo"
)

annot_cols <- c("Index", "NumberPSM", "Gene", "Protein", "Protein ID", "Protein Description")
sample_cols <- c(wt_pairs$igg, wt_pairs$ip, homo_pairs$igg, homo_pairs$ip)

dat <- df |>
  dplyr::select(dplyr::all_of(annot_cols), dplyr::all_of(sample_cols))

quant_mat <- dat |>
  dplyr::select(dplyr::all_of(sample_cols)) |>
  as.matrix()

quant_mat_imp <- impute_matrix_row_min(quant_mat)

dat_imp <- dplyr::bind_cols(
  dat |> dplyr::select(dplyr::all_of(annot_cols)),
  tibble::as_tibble(quant_mat_imp)
)

wt_ratio_tbl <- calculate_ip_igg_ratios(dat_imp, wt_pairs)
homo_ratio_tbl <- calculate_ip_igg_ratios(dat_imp, homo_pairs)

ratio_dat <- dplyr::bind_cols(
  dat_imp |> dplyr::select(dplyr::all_of(annot_cols)),
  wt_ratio_tbl,
  homo_ratio_tbl
)

wt_ratio_cols <- colnames(wt_ratio_tbl)
homo_ratio_cols <- colnames(homo_ratio_tbl)

final_interactome <- ratio_dat |>
  dplyr::rowwise() |>
  dplyr::mutate(
    WT_mean_log2_ratio = mean(c_across(dplyr::all_of(wt_ratio_cols)), na.rm = TRUE),
    WT_p_value = calculate_enrichment_pvalue(c_across(dplyr::all_of(wt_ratio_cols))),
    Homo_mean_log2_ratio = mean(c_across(dplyr::all_of(homo_ratio_cols)), na.rm = TRUE),
    Homo_p_value = calculate_enrichment_pvalue(c_across(dplyr::all_of(homo_ratio_cols)))
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    WT_q_value = p.adjust(WT_p_value, method = "BH"),
    Homo_q_value = p.adjust(Homo_p_value, method = "BH")
  ) |>
  classify_interactors(
    condition_a_col = "WT_mean_log2_ratio",
    condition_b_col = "Homo_mean_log2_ratio",
    condition_a_label = "WT-enriched",
    condition_b_label = "Homo-enriched"
  )

save_tsv(
  final_interactome,
  "results/BSN_interactome_combined_rowmin_imputed.tsv"
)
