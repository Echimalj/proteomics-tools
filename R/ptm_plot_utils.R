#' PTM plotting utilities
#'
#' @keywords internal
NULL

#' Build PTM lollipop plot
#'
#' @param ptm_df Processed PTM site table.
#' @param protein_len Protein length.
#' @param mutation_pos Optional mutation position.
#' @param title Plot title.
#'
#' @return ggplot object.
#' @export
plot_ptm_lollipop <- function(ptm_df,
                              protein_len,
                              mutation_pos = NULL,
                              mutation_label = NULL,
                              title = "PTM site map",
                              site_color = "#E74C3C") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required.", call. = FALSE)
  }

  ptm_long <- ptm_df |>
    dplyr::select(.data$site, .data$global_pos, .data$Peptide,
                  .data$`WT Intensity`, .data$`Homo Intensity`) |>
    tidyr::pivot_longer(
      cols = c(.data$`WT Intensity`, .data$`Homo Intensity`),
      names_to = "condition",
      values_to = "intensity"
    ) |>
    dplyr::mutate(
      condition = dplyr::recode(
        .data$condition,
        `WT Intensity` = "WT",
        `Homo Intensity` = "Homo"
      ),
      log_intensity = log10(.data$intensity + 1)
    )

  rng <- range(ptm_long$log_intensity, na.rm = TRUE)

  ptm_long <- ptm_long |>
    dplyr::mutate(
      y_base = 1,
      y_top = if (diff(rng) == 0) {
        1.8
      } else {
        1.2 + (.data$log_intensity - rng[1]) / diff(rng) * 1.2
      }
    )

  label_df <- ptm_long |>
    dplyr::group_by(.data$site, .data$global_pos, .data$condition) |>
    dplyr::summarise(
      y_top = max(.data$y_top, na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = data.frame(condition = unique(ptm_long$condition)),
      ggplot2::aes(xmin = 1, xmax = protein_len, ymin = 0, ymax = 1),
      fill = "grey90",
      color = "black"
    ) +
    ggplot2::geom_segment(
      data = ptm_long,
      ggplot2::aes(
        x = .data$global_pos,
        xend = .data$global_pos,
        y = .data$y_base,
        yend = .data$y_top
      ),
      color = site_color,
      alpha = 0.7
    ) +
    ggplot2::geom_point(
      data = ptm_long,
      ggplot2::aes(
        x = .data$global_pos,
        y = .data$y_top,
        size = .data$log_intensity
      ),
      color = site_color,
      alpha = 0.9
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(
        x = .data$global_pos,
        y = .data$y_top,
        label = .data$site
      ),
      size = 3,
      seed = 1
    ) +
    ggplot2::facet_grid(condition ~ .) +
    ggplot2::coord_cartesian(ylim = c(0, 3.1)) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      title = title,
      x = "Protein position",
      y = ""
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  if (!is.null(mutation_pos)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = mutation_pos,
        color = "#2ECC71",
        linetype = "dashed"
      )

    if (!is.null(mutation_label)) {
      p <- p +
        ggplot2::annotate(
          "text",
          x = mutation_pos,
          y = 2.95,
          label = mutation_label,
          color = "#27AE60",
          fontface = "bold",
          size = 3.5
        )
    }
  }

  p
}
