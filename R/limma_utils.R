#' limma utilities for proteomics
#'
#' @keywords internal
NULL

#' Build limma design matrix
#'
#' @param targets Sample metadata with group column.
#' @param group_col Group column.
#'
#' @return Design matrix.
#' @export
build_limma_design <- function(targets,
                               group_col = "group") {
  design <- stats::model.matrix(
    stats::as.formula(paste0("~ 0 + ", group_col)),
    data = targets
  )

  colnames(design) <- levels(targets[[group_col]])
  design
}

#' Fit limma contrasts
#'
#' @param expr Numeric matrix.
#' @param design Design matrix.
#' @param contrast_matrix limma contrast matrix.
#'
#' @return eBayes fit object.
#' @export
fit_limma_contrasts <- function(expr,
                                design,
                                contrast_matrix) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required.", call. = FALSE)
  }

  fit <- limma::lmFit(expr, design)
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  limma::eBayes(fit2)
}

#' Extract limma result table
#'
#' @param fit limma fit object.
#' @param coef Contrast name.
#'
#' @return A tibble/data frame.
#' @export
extract_limma_results <- function(fit, coef) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required.", call. = FALSE)
  }

  limma::topTable(
    fit,
    coef = coef,
    number = Inf,
    sort.by = "P"
  ) |>
    tibble::rownames_to_column("row_id")
}
