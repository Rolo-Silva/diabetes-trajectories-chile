# R/02_trajectories.R
# Trajectory plots for manuscript outputs:
#   - SF3: national mean trends (EGCC + DRSC, two panels)
#   - Figure 1B: EGCC class-specific trajectories
#   - Figure 2B: DRSC class-specific trajectories


# -------------------------------------------------------------------------
# make_mean_trends_plot()
# Two-panel plot of national mean coverage trends (EGCC + DRSC), 2011-2023,
# with 95% CI error bars and 80% WHO benchmark line.
# -------------------------------------------------------------------------

make_mean_trends_plot <- function(coverage_data,
                                  out_path = "outputs/figures/sf3_mean_trends.png",
                                  width = 9, height = 5, dpi = 300) {
  
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  summary_data <- coverage_data |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      mean_egcc = mean(dgcc, na.rm = TRUE),
      se_egcc   = sd(dgcc, na.rm = TRUE) / sqrt(sum(!is.na(dgcc))),
      mean_drsc = mean(drsc, na.rm = TRUE),
      se_drsc   = sd(drsc, na.rm = TRUE) / sqrt(sum(!is.na(drsc))),
      .groups   = "drop"
    ) |>
    dplyr::mutate(
      lower_egcc = mean_egcc - 1.96 * se_egcc,
      upper_egcc = mean_egcc + 1.96 * se_egcc,
      lower_drsc = mean_drsc - 1.96 * se_drsc,
      upper_drsc = mean_drsc + 1.96 * se_drsc
    )
  
  plot_panel <- function(data, mean_var, lower_var, upper_var, y_label) {
    ggplot2::ggplot(data, ggplot2::aes(x = year, y = .data[[mean_var]])) +
      ggplot2::geom_line(linewidth = 1, colour = "red") +
      ggplot2::geom_point(size = 1.5, colour = "red") +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[[lower_var]], ymax = .data[[upper_var]]),
        width = 0.25, colour = "red"
      ) +
      ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", colour = "gray50") +
      ggplot2::scale_x_continuous(
        breaks = seq(0, 12, by = 2),
        labels = seq(2011, 2023, by = 2)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2),
        labels = scales::percent
      ) +
      ggplot2::xlab("Year") +
      ggplot2::ylab(y_label) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text        = ggplot2::element_text(size = 10),
        axis.title       = ggplot2::element_text(size = 10, face = "bold"),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
  }
  
  panel_a <- plot_panel(summary_data, "mean_egcc", "lower_egcc", "upper_egcc",
                        "Effective glycaemic control coverage")
  panel_b <- plot_panel(summary_data, "mean_drsc", "lower_drsc", "upper_drsc",
                        "Diabetic retinopathy screening coverage")
  
  combined <- panel_a | panel_b
  
  ggplot2::ggsave(out_path, plot = combined,
                  width = width, height = height, dpi = dpi, units = "in")
  
  return(out_path)
}


# -------------------------------------------------------------------------
# make_class_trajectory_summary()
# Internal helper: aggregates coverage data by year and class, returning
# mean, SD, n, and 95% CI per cell.
# -------------------------------------------------------------------------

make_class_trajectory_summary <- function(coverage_data_with_classes, outcome, class_var) {
  coverage_data_with_classes |>
    dplyr::group_by(year, {{ class_var }}) |>
    dplyr::summarise(
      mean_coverage = mean({{ outcome }}, na.rm = TRUE),
      sd_coverage   = sd({{ outcome }},   na.rm = TRUE),
      n_coverage    = sum(!is.na({{ outcome }})),
      se_mean       = sd_coverage / sqrt(n_coverage),
      lower_ci      = mean_coverage - 1.96 * se_mean,
      upper_ci      = mean_coverage + 1.96 * se_mean,
      .groups       = "drop"
    ) |>
    dplyr::rename(class = {{ class_var }})
}


# -------------------------------------------------------------------------
# make_class_trajectory_plot()
# Internal helper: builds a single class-specific trajectory plot from
# pre-computed summary data.
# -------------------------------------------------------------------------

make_class_trajectory_plot <- function(summary_data, class_labels, y_axis_label,
                                       panel_letter = NULL) {
  
  n_classes <- length(class_labels)
  
  p <- summary_data |>
    dplyr::mutate(class = factor(class, levels = seq_len(n_classes),
                                 labels = class_labels)) |>
    ggplot2::ggplot(ggplot2::aes(x = year, y = mean_coverage, colour = class)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower_ci, ymax = upper_ci), width = 0.25) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", colour = "gray50") +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 12, by = 2),
      labels = seq(2011, 2023, by = 2)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = scales::percent
    ) +
    ggplot2::xlab("Year") +
    ggplot2::ylab(y_axis_label) +
    ggplot2::labs(colour = "Latent Class") +
    ggsci::scale_color_lancet() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title            = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.title          = ggplot2::element_text(size = 12, face = "bold"),
      axis.text             = ggplot2::element_text(size = 10),
      axis.title            = ggplot2::element_text(size = 10, face = "bold"),
      legend.position       = c(0.5, 0.98),
      legend.justification  = c(0.5, 1),
      legend.direction      = "horizontal",
      legend.background     = ggplot2::element_rect(fill = scales::alpha("white", 0.7), colour = NA),
      legend.box.background = ggplot2::element_blank(),
      legend.margin         = ggplot2::margin(0, 0, 0, 0),
      legend.key.height     = grid::unit(8, "pt"),
      legend.key.width      = grid::unit(16, "pt"),
      panel.grid.major      = ggplot2::element_blank(),
      panel.grid.minor      = ggplot2::element_blank(),
      legend.text           = ggplot2::element_text(size = 10)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE,
                                     title.position = "top", title.hjust = 0)
    )
  
  if (!is.null(panel_letter)) {
    p <- p + ggplot2::annotate(
      "text", x = max(summary_data$year), y = 1, label = panel_letter,
      hjust = 0.8, vjust = 1.2, size = 6, fontface = "bold"
    )
  }
  
  p
}


# -------------------------------------------------------------------------
# make_class_trajectories_plot()
# Main entry point: builds summary, plot, saves PNG. Returns file path
# for targets file tracking.
# -------------------------------------------------------------------------

make_class_trajectories_plot <- function(coverage_data_with_classes,
                                         outcome, class_var,
                                         class_labels, y_axis_label,
                                         out_path,
                                         panel_letter = NULL,
                                         width = 7, height = 5, dpi = 300) {
  
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  summary_data <- make_class_trajectory_summary(
    coverage_data_with_classes,
    outcome   = {{ outcome }},
    class_var = {{ class_var }}
  )
  
  p <- make_class_trajectory_plot(
    summary_data = summary_data,
    class_labels = class_labels,
    y_axis_label = y_axis_label,
    panel_letter = panel_letter
  )
  
  ggplot2::ggsave(out_path, plot = p,
                  width = width, height = height, dpi = dpi, units = "in")
  
  return(out_path)
}