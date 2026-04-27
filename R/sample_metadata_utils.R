#' Sample metadata utilities
#'
#' @keywords internal
NULL

#' Build proteomics sample metadata
#'
#' @param sample_cols Character vector of sample column names.
#' @param groups Character vector of group labels.
#' @param group_levels Optional group order.
#'
#' @return A tibble with sample and group columns.
#' @export
build_sample_metadata <- function(sample_cols,
                                  groups,
                                  group_levels = NULL) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.", call. = FALSE)
  }

  if (length(sample_cols) != length(groups)) {
    stop("sample_cols and groups must have the same length.", call. = FALSE)
  }

  if (is.null(group_levels)) {
    group_levels <- unique(groups)
  }

  tibble::tibble(
    sample = sample_cols,
    group = factor(groups, levels = group_levels)
  )
}

#' Build paired IP/IgG metadata table
#'
#' @param replicate Replicate IDs.
#' @param igg IgG sample columns.
#' @param ip IP sample columns.
#' @param condition Condition label.
#'
#' @return A tibble.
#' @export
build_ip_pairs <- function(replicate,
                           igg,
                           ip,
                           condition) {
  tibble::tibble(
    replicate = replicate,
    condition = condition,
    igg = igg,
    ip = ip
  )
}
