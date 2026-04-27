#' Export utilities
#'
#' @keywords internal
NULL

#' Save table as CSV
#'
#' @param x Data frame.
#' @param file Output file.
#'
#' @return Invisibly returns file path.
#' @export
save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file)
  invisible(file)
}

#' Save table as TSV
#'
#' @param x Data frame.
#' @param file Output file.
#'
#' @return Invisibly returns file path.
#' @export
save_tsv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, file)
  invisible(file)
}

#' Save ggplot
#'
#' @param p ggplot object.
#' @param file Output file.
#' @param width Width.
#' @param height Height.
#'
#' @return Invisibly returns file path.
#' @export
save_plot <- function(p, file, width = 8, height = 6) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file, p, width = width, height = height)
  invisible(file)
}
