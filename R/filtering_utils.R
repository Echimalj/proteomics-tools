#' Filtering utilities for proteomics matrices
#'
#' @keywords internal
NULL

#' Filter proteins by completeness across paired groups
#'
#' @param expr Numeric matrix.
#' @param group_cols Named list of sample column vectors.
#' @param min_observed Minimum observed values per group.
#'
#' @return Logical vector.
#' @export
filter_by_group_completeness <- function(expr,
                                         group_cols,
                                         min_observed = 2) {
  keep <- rep(TRUE, nrow(expr))

  for (nm in names(group_cols)) {
    cols <- group_cols[[nm]]
    missing_cols <- setdiff(cols, colnames(expr))

    if (length(missing_cols) > 0) {
      stop(
        "Missing columns for group ", nm, ": ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    keep <- keep & rowSums(!is.na(expr[, cols, drop = FALSE])) >= min_observed
  }

  keep
}

#' Summarize filtering
#'
#' @param before_n Number before filtering.
#' @param after_n Number after filtering.
#'
#' @return A data frame.
#' @export
summarize_filtering <- function(before_n, after_n) {
  data.frame(
    before = before_n,
    after = after_n,
    removed = before_n - after_n,
    retained_fraction = after_n / before_n
  )
}
