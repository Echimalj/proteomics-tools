source("R/load_fragpipe_utils.R")
source("R/ptm_filter_utils.R")
source("R/ptm_plot_utils.R")
source("R/export_utils.R")

library(tidyverse)

wt_psm <- read_fragpipe_psm_table("data/psm_WT.tsv", condition = "WT")
homo_psm <- read_fragpipe_psm_table("data/psm_Homo.tsv", condition = "Homo")

psm_all <- dplyr::bind_rows(wt_psm, homo_psm)

bsn_psm <- filter_modified_target_psms(
  psm_df = psm_all,
  gene = "Bsn",
  qvalue_cutoff = 0.01
)

ptm_presence <- summarize_ptm_presence(bsn_psm) |>
  classify_ptm_type(ptm_col = "ptm_id")

save_tsv(ptm_presence, "results/BSN_PTM_presence_all.tsv")

phospho <- read_fragpipe_site_table("data/combined_site_STY_79.9663.tsv")

phospho_bsn <- process_target_site_table(
  site_df = phospho,
  gene = "Bsn",
  site_regex = "S\\d+"
)

save_tsv(phospho_bsn, "results/BSN_phosphorylation_sites.tsv")

p_phospho <- plot_ptm_lollipop(
  ptm_df = phospho_bsn,
  protein_len = 3942,
  mutation_pos = 3882,
  mutation_label = "P3882A Mutation",
  title = "Bassoon phosphorylation map: WT vs Homo"
)

save_plot(p_phospho, "figures/BSN_phosphorylation_lollipop.pdf", width = 12, height = 10)

acetyl <- read_fragpipe_site_table("data/combined_site_Kn_42.0106.tsv")

acetyl_bsn <- process_target_site_table(
  site_df = acetyl,
  gene = "Bsn",
  site_regex = "[SK]\\d+"
)

save_tsv(acetyl_bsn, "results/BSN_acetylation_sites.tsv")
