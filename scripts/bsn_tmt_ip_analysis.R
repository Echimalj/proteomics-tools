source("R/load_fragpipe_utils.R")
source("R/sample_metadata_utils.R")
source("R/filtering_utils.R")
source("R/limma_utils.R")
source("R/interactome_utils.R")
source("R/plotting_utils.R")
source("R/export_utils.R")

library(tidyverse)
library(limma)

df <- read_fragpipe_protein_table(
  "data/abundance_protein_MD.tsv"
)

sample_cols <- c(
  "TP652-WT-IgG", "BN141-WT-IgG", "BN368-WT-IgG",
  "TP652-WT-BSN_IP", "BN141-WT-BSN_IP", "BN368-WT-BSN_IP",
  "BN207-BSNKIHomo-IgG", "BN206-BSNKIHomo-IgG", "BN214-BSNKIHomo-IgG",
  "BN207-BSNKIHomo-BSNIP", "BN206-BSNKIHomo-BSNIP", "BN214-BSNKIHomo-BSNIP"
)

dat <- df |>
  dplyr::select(Gene, `Protein Description`, `Protein ID`, dplyr::all_of(sample_cols))

expr <- dat |>
  dplyr::select(dplyr::all_of(sample_cols)) |>
  as.matrix()

rownames(expr) <- make.unique(paste(dat$Gene, dat$`Protein ID`, sep = "|"))

targets <- build_sample_metadata(
  sample_cols = sample_cols,
  groups = c(
    "WT_IgG", "WT_IgG", "WT_IgG",
    "WT_IP", "WT_IP", "WT_IP",
    "Homo_IgG", "Homo_IgG", "Homo_IgG",
    "Homo_IP", "Homo_IP", "Homo_IP"
  ),
  group_levels = c("WT_IgG", "WT_IP", "Homo_IgG", "Homo_IP")
)

keep <- filter_by_group_completeness(
  expr,
  group_cols = list(
    WT_IgG = c("TP652-WT-IgG", "BN141-WT-IgG", "BN368-WT-IgG"),
    WT_IP = c("TP652-WT-BSN_IP", "BN141-WT-BSN_IP", "BN368-WT-BSN_IP"),
    Homo_IgG = c("BN207-BSNKIHomo-IgG", "BN206-BSNKIHomo-IgG", "BN214-BSNKIHomo-IgG"),
    Homo_IP = c("BN207-BSNKIHomo-BSNIP", "BN206-BSNKIHomo-BSNIP", "BN214-BSNKIHomo-BSNIP")
  ),
  min_observed = 2
)

expr_filt <- expr[keep, , drop = FALSE]
dat_filt <- dat[keep, , drop = FALSE]

design <- build_limma_design(targets)

contrast_matrix <- limma::makeContrasts(
  WT_IP_vs_WT_IgG = WT_IP - WT_IgG,
  Homo_IP_vs_Homo_IgG = Homo_IP - Homo_IgG,
  Delta_Interaction = (WT_IP - WT_IgG) - (Homo_IP - Homo_IgG),
  WT_IP_vs_Homo_IP = WT_IP - Homo_IP,
  levels = design
)

fit <- fit_limma_contrasts(expr_filt, design, contrast_matrix)

res_main <- extract_limma_results(fit, coef = "Delta_Interaction")

annot_tbl <- dat_filt |>
  dplyr::mutate(row_id = rownames(expr_filt)) |>
  dplyr::select(row_id, Gene, `Protein ID`, `Protein Description`)

summary_tbl <- dat_filt |>
  dplyr::mutate(
    WT_IgG_mean = rowMeans(dplyr::select(., `TP652-WT-IgG`, `BN141-WT-IgG`, `BN368-WT-IgG`), na.rm = TRUE),
    WT_IP_mean = rowMeans(dplyr::select(., `TP652-WT-BSN_IP`, `BN141-WT-BSN_IP`, `BN368-WT-BSN_IP`), na.rm = TRUE),
    Homo_IgG_mean = rowMeans(dplyr::select(., `BN207-BSNKIHomo-IgG`, `BN206-BSNKIHomo-IgG`, `BN214-BSNKIHomo-IgG`), na.rm = TRUE),
    Homo_IP_mean = rowMeans(dplyr::select(., `BN207-BSNKIHomo-BSNIP`, `BN206-BSNKIHomo-BSNIP`, `BN214-BSNKIHomo-BSNIP`), na.rm = TRUE),
    WT_enrichment = WT_IP_mean - WT_IgG_mean,
    Homo_enrichment = Homo_IP_mean - Homo_IgG_mean,
    delta_enrichment = WT_enrichment - Homo_enrichment,
    row_id = rownames(expr_filt)
  ) |>
  dplyr::select(row_id, WT_enrichment, Homo_enrichment, delta_enrichment)

res_main <- res_main |>
  dplyr::left_join(annot_tbl, by = "row_id") |>
  dplyr::left_join(summary_tbl, by = "row_id") |>
  dplyr::relocate(Gene, `Protein ID`, `Protein Description`)

p <- plot_volcano(res_main, x_col = "logFC", p_col = "P.Value", label_col = "Gene")

save_csv(res_main, "results/BSN_IP_limma_results_full.csv")
save_plot(p, "figures/BSN_IP_limma_volcano.pdf", width = 7, height = 6)
