#' Missing-value imputation utilities
#'
#' @keywords internal
NULL

#' Row-min imputation
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector.
#' @export
impute_row_min <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }

  row_min <- min(x, na.rm = TRUE)
  x[is.na(x)] <- row_min
  x
}

#' Apply row-min imputation to matrix
#'
#' @param mat Numeric matrix.
#'
#' @return Numeric matrix.
#' @export
impute_matrix_row_min <- function(mat) {
  out <- t(apply(mat, 1, impute_row_min))
  colnames(out) <- colnames(mat)
  rownames(out) <- rownames(mat)
  out
}
