# R/06_cross_class.R
# -------------------------------------------------------------------------
# Cross-classification analysis of EGCC × DRSC trajectory classes.
#
# Modal class assignment is derived independently for each outcome from
# posterior class probabilities (prob1..prob4) via max.col(), replicating
# the inferential path used in the original analysis. This is
# mathematically equivalent to lcmm's $pprob$class output but makes the
# assignment rule explicit and defensible under peer review.
#
# Computes:
#   - 4×4 cross-tabulation with counts and row percentages
#   - Pearson chi-square (asymptotic + Monte Carlo p-value, B=10000)
#   - Standardised residuals
#   - Cramér's V
#   - Best-match summaries (row-wise and column-wise)
#   - Heatmap figure (Supplementary Figure 9)
#
# Returns a named list with all components plus the file path of the
# saved heatmap.
# -------------------------------------------------------------------------

make_cross_classification <- function(coverage_demographics_class,
                                      out_path  = "outputs/figures/sf9_cross_classification.png",
                                      width = 9, height = 6, dpi = 300,
                                      chi_sim_B = 10000) {
  
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  # ---- Labels for the 4 classes of each outcome ----
  
  egcc_labels <- c("Lowest", "Highest", "Stable upper medium", "Stable medium")
  drsc_labels <- c("Increasing", "Stable medium",
                   "Predominantly highest", "Lowest")
  
  prob_cols_dgcc <- c("prob1_dgcc", "prob2_dgcc", "prob3_dgcc", "prob4_dgcc")
  prob_cols_drsc <- c("prob1_drsc", "prob2_drsc", "prob3_drsc", "prob4_drsc")
  
  stopifnot(all(prob_cols_dgcc %in% names(coverage_demographics_class)))
  stopifnot(all(prob_cols_drsc %in% names(coverage_demographics_class)))
  
  # ---- Version B: derive modal class from posterior probabilities ----
  
  probs_dgcc <- as.matrix(coverage_demographics_class[, prob_cols_dgcc])
  probs_drsc <- as.matrix(coverage_demographics_class[, prob_cols_drsc])
  
  classes_df <- coverage_demographics_class |>
    dplyr::select(id) |>
    dplyr::mutate(
      egcc_class_num   = max.col(probs_dgcc, ties.method = "first"),
      drsc_class_num   = max.col(probs_drsc, ties.method = "first"),
      egcc_class_label = factor(egcc_class_num, levels = 1:4, labels = egcc_labels),
      drsc_class_label = factor(drsc_class_num, levels = 1:4, labels = drsc_labels)
    ) |>
    tidyr::drop_na(egcc_class_label, drsc_class_label)
  
  # ---- 4×4 cross-tabulation ----
  
  tab     <- table(classes_df$egcc_class_label, classes_df$drsc_class_label)
  row_pct <- prop.table(tab, margin = 1) * 100
  col_pct <- prop.table(tab, margin = 2) * 100
  
  # ---- Final table: n (row%) with row totals + column total ----
  
  tab_n_pct <- matrix(
    paste0(tab, " (", sprintf("%.1f", row_pct), "%)"),
    nrow = nrow(tab),
    ncol = ncol(tab),
    dimnames = dimnames(tab)
  )
  
  tab_final <- as.data.frame.matrix(tab_n_pct)
  tab_final$Row_Total_n <- rowSums(tab)
  
  total_row <- c(as.character(colSums(tab)), sum(tab))
  names(total_row) <- colnames(tab_final)
  tab_final <- rbind(tab_final, Total_n = total_row)
  
  # ---- Chi-square: asymptotic + Monte Carlo ----
  
  chisq_asym <- suppressWarnings(stats::chisq.test(tab))
  chisq_mc   <- suppressWarnings(
    stats::chisq.test(tab, simulate.p.value = TRUE, B = chi_sim_B)
  )
  
  # ---- Cramér's V ----
  
  n <- sum(tab)
  k <- min(nrow(tab) - 1, ncol(tab) - 1)
  cramers_v <- sqrt(as.numeric(chisq_asym$statistic) / (n * k))
  
  # ---- Best-match summaries ----
  
  best_match_by_egcc <- data.frame(
    egcc_class      = rownames(row_pct),
    best_drsc_class = apply(row_pct, 1, function(x) names(which.max(x))),
    percentage      = round(apply(row_pct, 1, max), 1),
    stringsAsFactors = FALSE
  )
  
  best_match_by_drsc <- data.frame(
    drsc_class      = colnames(col_pct),
    best_egcc_class = apply(col_pct, 2, function(x) names(which.max(x))),
    percentage      = round(apply(col_pct, 2, max), 1),
    stringsAsFactors = FALSE
  )
  
  # ---- Heatmap data ----
  
  heatmap_df <- as.data.frame(as.table(tab))
  colnames(heatmap_df) <- c("EGCC", "DRSC", "n")
  heatmap_df$row_pct <- as.vector(row_pct)
  heatmap_df$label   <- paste0(
    heatmap_df$n, "\n(",
    sprintf("%.1f", heatmap_df$row_pct), "%)"
  )
  
  heatmap_df$EGCC <- factor(heatmap_df$EGCC, levels = rownames(tab))
  heatmap_df$DRSC <- factor(heatmap_df$DRSC, levels = colnames(tab))
  
  # ---- Heatmap plot ----
  
  p <- ggplot2::ggplot(heatmap_df,
                       ggplot2::aes(x = DRSC, y = EGCC, fill = row_pct)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.8) +
    ggplot2::scale_fill_gradient(low = "grey95", high = "steelblue",
                                 name = "Row %") +
    ggplot2::labs(
      title    = "Cross-classification of EGCC and DRSC trajectories",
      subtitle = "Cells show number of municipalities and row percentage",
      x = "DRSC trajectory class",
      y = "EGCC trajectory class"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid      = ggplot2::element_blank(),
      axis.text.x     = ggplot2::element_text(angle = 25, hjust = 1),
      plot.title      = ggplot2::element_text(face = "bold")
    )
  
  ggplot2::ggsave(out_path, plot = p,
                  width = width, height = height,
                  dpi = dpi, units = "in")
  
  # ---- Return all components ----
  
  list(
    classes_df         = classes_df,
    tab                = tab,
    tab_with_pct       = tab_final,
    row_percentages    = row_pct,
    col_percentages    = col_pct,
    expected_counts    = chisq_asym$expected,
    std_residuals      = chisq_asym$stdres,
    chisq_asymptotic   = chisq_asym,
    chisq_monte_carlo  = chisq_mc,
    cramers_v          = cramers_v,
    best_match_by_egcc = best_match_by_egcc,
    best_match_by_drsc = best_match_by_drsc,
    plot_path          = out_path
  )
}