#' Load FragPipe output tables
#'
#' @keywords internal
NULL

#' Read FragPipe protein abundance table
#'
#' @param file Path to FragPipe protein abundance table.
#'
#' @return A tibble.
#' @export
read_fragpipe_protein_table <- function(file) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required.", call. = FALSE)
  }

  readr::read_tsv(file, show_col_types = FALSE)
}

#' Read FragPipe PSM table
#'
#' @param file Path to psm.tsv.
#' @param condition Optional condition label to add.
#'
#' @return A tibble.
#' @export
read_fragpipe_psm_table <- function(file, condition = NULL) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required.", call. = FALSE)
  }

  x <- readr::read_tsv(file, show_col_types = FALSE)

  if (!is.null(condition)) {
    x$condition <- condition
  }

  x
}

#' Read FragPipe PTM site table
#'
#' @param file Path to combined_site table.
#'
#' @return A tibble.
#' @export
read_fragpipe_site_table <- function(file) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required.", call. = FALSE)
  }

  readr::read_tsv(file, show_col_types = FALSE)
}
