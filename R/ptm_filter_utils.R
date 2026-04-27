#' PTM filtering utilities
#'
#' @keywords internal
NULL

#' Filter PSM table for modified target protein peptides
#'
#' @param psm_df PSM table.
#' @param gene Target gene.
#' @param qvalue_cutoff Q-value cutoff.
#'
#' @return Filtered PSM table.
#' @export
filter_modified_target_psms <- function(psm_df,
                                        gene,
                                        qvalue_cutoff = 0.01) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required.", call. = FALSE)
  }

  psm_df |>
    dplyr::filter(
      .data$Gene == gene,
      .data$Qvalue <= qvalue_cutoff,
      .data$`Is Decoy` == FALSE,
      .data$`Is Contaminant` == FALSE,
      !is.na(.data$`Assigned Modifications`),
      .data$`Assigned Modifications` != ""
    ) |>
    dplyr::mutate(
      run_id = stringr::str_remove(.data$Spectrum, "\\.\\d+\\.\\d+\\.\\d+$"),
      ptm_id = paste(
        .data$`Modified Peptide`,
        .data$`Assigned Modifications`,
        .data$`Best Positions`,
        sep = " | "
      )
    )
}

#' Summarize PTM presence by condition
#'
#' @param psm_df Filtered PSM table.
#'
#' @return PTM summary table.
#' @export
summarize_ptm_presence <- function(psm_df) {
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }

  psm_df |>
    dplyr::group_by(.data$ptm_id, .data$condition) |>
    dplyr::summarise(
      n_psm = dplyr::n(),
      n_runs = dplyr::n_distinct(.data$run_id),
      mean_intensity = mean(.data$Intensity, na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = .data$condition,
      values_from = c(.data$n_psm, .data$n_runs, .data$mean_intensity),
      values_fill = 0
    )
}

#' Classify PTM type
#'
#' @param ptm_df PTM table.
#' @param ptm_col Column containing PTM annotation.
#'
#' @return PTM table with PTM_type.
#' @export
classify_ptm_type <- function(ptm_df,
                              ptm_col = "ptm_id") {
  ptm_df |>
    dplyr::mutate(
      PTM_type = dplyr::case_when(
        stringr::str_detect(.data[[ptm_col]], "57\\.0214") ~ "Carbamidomethyl_C",
        stringr::str_detect(.data[[ptm_col]], "15\\.9949") ~ "Oxidation_M",
        stringr::str_detect(.data[[ptm_col]], "42\\.0106") ~ "Nterm_Acetyl",
        TRUE ~ "Other"
      )
    )
}

#' Process target protein PTM site table
#'
#' @param site_df FragPipe combined site table.
#' @param gene Target gene.
#' @param site_regex Regex for site extraction.
#'
#' @return Processed site table.
#' @export
process_target_site_table <- function(site_df,
                                      gene,
                                      site_regex = "S\\d+") {
  site_df |>
    dplyr::filter(.data$Gene == gene) |>
    dplyr::mutate(
      site = stringr::str_extract(.data$Index, site_regex),
      global_pos = as.numeric(stringr::str_extract(.data$site, "\\d+")),
      log2FC = log2((.data$`WT Intensity` + 1) / (.data$`Homo Intensity` + 1)),
      Class = dplyr::case_when(
        .data$log2FC > 1 ~ "WT-enriched",
        .data$log2FC < -1 ~ "Homo-enriched",
        TRUE ~ "Shared"
      )
    )
}
