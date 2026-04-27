#' Interactome analysis utilities
#'
#' @keywords internal
NULL

#' Calculate IP minus IgG log2 ratios
#'
#' @param dat Data frame with abundance columns.
#' @param pair_tbl Table with replicate, igg, and ip columns.
#' @param suffix Output suffix.
#'
#' @return A tibble of paired ratios.
#' @export
calculate_ip_igg_ratios <- function(dat,
                                    pair_tbl,
                                    suffix = "log2ratio") {
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package 'purrr' is required.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }

  out <- purrr::map2_dfc(
    pair_tbl$ip,
    pair_tbl$igg,
    ~ dat[[.x]] - dat[[.y]]
  )

  colnames(out) <- paste0(
    pair_tbl$replicate,
    "_",
    pair_tbl$condition,
    "_",
    suffix
  )

  out
}

#' One-sided enrichment p-value
#'
#' @param x Numeric vector.
#'
#' @return p-value.
#' @export
calculate_enrichment_pvalue <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) < 2) {
    return(NA_real_)
  }

  tryCatch(
    stats::t.test(x, alternative = "greater", mu = 0)$p.value,
    error = function(e) NA_real_
  )
}

#' Classify condition-enriched interactors
#'
#' @param df Data frame.
#' @param condition_a_col Mean enrichment column for condition A.
#' @param condition_b_col Mean enrichment column for condition B.
#' @param min_enrichment Minimum enrichment.
#' @param min_delta Minimum difference.
#'
#' @return Data frame with Class column.
#' @export
classify_interactors <- function(df,
                                 condition_a_col,
                                 condition_b_col,
                                 condition_a_label = "ConditionA-enriched",
                                 condition_b_label = "ConditionB-enriched",
                                 min_enrichment = 0.5,
                                 min_delta = 0.5) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }

  a <- df[[condition_a_col]]
  b <- df[[condition_b_col]]
  delta <- a - b

  df |>
    dplyr::mutate(
      delta_enrichment = delta,
      Class = dplyr::case_when(
        a > min_enrichment & delta > min_delta ~ condition_a_label,
        b > min_enrichment & delta < -min_delta ~ condition_b_label,
        TRUE ~ "Shared/weak"
      )
    )
}
