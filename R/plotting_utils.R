#' Plotting utilities
#'
#' @keywords internal
NULL

#' Basic volcano plot
#'
#' @param df Result table.
#' @param x_col X-axis column.
#' @param p_col P-value column.
#' @param label_col Label column.
#' @param label_top_n Number of labels.
#'
#' @return ggplot object.
#' @export
plot_volcano <- function(df,
                         x_col = "logFC",
                         p_col = "P.Value",
                         label_col = "Gene",
                         label_top_n = 10) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required.", call. = FALSE)
  }

  df$neg_log10_p <- -log10(df[[p_col]])

  top <- df |>
    dplyr::arrange(.data[[p_col]]) |>
    head(label_top_n)

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_col]], y = .data$neg_log10_p)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggrepel::geom_text_repel(
      data = top,
      ggplot2::aes(label = .data[[label_col]])
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = x_col,
      y = "-log10(p-value)"
    )
}
