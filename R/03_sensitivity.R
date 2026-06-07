# R/03_sensitivity.R
# -------------------------------------------------------------------------
# Sensitivity analyses for the retained 4-class cubic NRE trajectory models.
#
# Two independent validation stages:
#   - Stage 1: municipality subsampling validation (70%, 80%, 90% of IDs).
#     Tests how trajectory shapes change when entire municipalities are
#     randomly excluded.
#
#   - Stage 2: municipality-year observation perturbation (5%, 10%, 15%, 20%).
#     Tests how trajectory shapes change when individual time-point
#     observations are randomly removed.
#
# Both stages REFIT the selected model on each perturbed dataset and
# compare predicted mean trajectories to the full-sample reference via
# RMSE and max absolute deviation, after class re-alignment.
#
# DEFAULTS: 1 seed (fast pipeline). Override seed_vec via _targets.R
# parameter to run 20-seed production validation.
# -------------------------------------------------------------------------


# =========================================================================
# Shared helpers (all prefixed .sens_*)
# =========================================================================

# Reference model selector for the saved full-sample fits
.sens_ref_model_for <- function(var_name, all_named_models) {
  if (var_name == "dgcc") {
    all_named_models[["4class_cubic_nre_dgcc_model"]]
  } else if (var_name == "drsc") {
    all_named_models[["4class_cubic_nre_drsc_model"]]
  } else {
    NULL
  }
}

# Find the time variable in a data frame
.sens_find_time_var <- function(df,
                                prefer = c("year", "ano", "time", "t", "Time", "timevar")) {
  v <- intersect(prefer, names(df))
  if (length(v)) v[1] else
    stop("Time variable not found in data (tried: ",
         paste(prefer, collapse = ", "), ").")
}

# Safe numeric coercion of time
.sens_as_numeric_time <- function(x) {
  if (inherits(x, "Date"))   return(as.numeric(x))
  if (inherits(x, "POSIXt")) return(as.numeric(x))
  suppressWarnings({
    xn <- as.numeric(x)
    if (all(is.finite(xn) | is.na(xn))) return(xn)
  })
  stop("time_var must be numeric/Date/POSIX; couldn't coerce safely.")
}

# Build polynomial terms RHS
.sens_poly_terms <- function(time_var = "year", deg = 3) {
  deg <- as.integer(deg)
  stopifnot(deg %in% 1:3)
  paste(
    c(time_var,
      if (deg >= 2) sprintf("I(%s^2)", time_var),
      if (deg >= 3) sprintf("I(%s^3)", time_var)),
    collapse = " + "
  )
}

# All permutations of 1:K (no extra package required)
.sens_perms <- function(v) {
  if (length(v) == 1) return(list(v))
  out <- list()
  for (i in seq_along(v)) {
    rest <- v[-i]
    for (p in .sens_perms(rest)) out[[length(out) + 1]] <- c(v[i], p)
  }
  out
}

# Build mean predicted trajectory curve (class_rank x yearv x pred)
# from a fitted lcmm model and a corresponding data frame.
.sens_pred_curve_from_model <- function(model, data_full, time_var = "year", K = 4) {
  
  pred <- model$pred
  if (is.null(pred)) return(NULL)
  
  pred_cols <- grep("^pred_m\\d+$", names(pred), value = TRUE)
  if (!length(pred_cols)) return(NULL)
  
  id_col <- if ("id" %in% names(pred)) "id" else names(pred)[1]
  
  data_aligned <- data_full |>
    dplyr::arrange(.data[[id_col]], .data[[time_var]]) |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::mutate(row_id = dplyr::row_number()) |>
    dplyr::ungroup()
  
  pred_aligned <- pred |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::mutate(row_id = dplyr::row_number()) |>
    dplyr::ungroup()
  
  combined <- dplyr::left_join(
    data_aligned |>
      dplyr::select(dplyr::all_of(id_col), row_id, dplyr::all_of(time_var)),
    pred_aligned |>
      dplyr::select(dplyr::all_of(id_col), row_id, dplyr::all_of(pred_cols)),
    by = c(id_col, "row_id")
  )
  
  combined |>
    tidyr::pivot_longer(
      dplyr::all_of(pred_cols),
      names_to = "class_label",
      values_to = "pred"
    ) |>
    dplyr::mutate(
      class = readr::parse_number(class_label),
      pred  = as.numeric(pred),
      yearv = .data[[time_var]]
    ) |>
    dplyr::group_by(class, yearv) |>
    dplyr::summarise(pred = mean(pred, na.rm = TRUE), .groups = "drop") |>
    dplyr::rename(class_rank = class)
}

# Align a perturbed run's class_rank to the reference by minimising
# total RMSE over class-specific mean curves.
.sens_align_to_ref_rmse <- function(run_curve, ref_curve, K = 4) {
  
  classes   <- seq_len(K)
  best_perm <- NULL
  best_val  <- Inf
  
  for (perm in .sens_perms(classes)) {
    
    total <- 0
    ok    <- TRUE
    
    for (r in classes) {
      
      j <- dplyr::inner_join(
        dplyr::filter(run_curve, class_rank == perm[r]),
        dplyr::filter(ref_curve, class_rank == r),
        by = "yearv",
        suffix = c("", "_ref")
      )
      
      if (nrow(j) == 0) {
        ok <- FALSE
        break
      }
      
      total <- total + sqrt(mean((j$pred - j$pred_ref)^2, na.rm = TRUE))
    }
    
    if (ok && total < best_val) {
      best_val  <- total
      best_perm <- perm
    }
  }
  
  inv <- integer(K)
  for (r in classes) inv[best_perm[r]] <- r
  
  dplyr::mutate(run_curve, class_rank = inv[class_rank])
}

# Refit the selected 4-class cubic NRE model on a supplied data frame.
# Returns list with model + pred_curve + fit_ok flag.
.sens_fit_selected_model <- function(var_name, data_in, K = 4L,
                                     time_var = "year",
                                     poly_degree = 3,
                                     seed = 123) {
  
  stopifnot(var_name %in% c("dgcc", "drsc"))
  
  df0 <- data_in |>
    tidyr::drop_na(id, dplyr::all_of(var_name), dplyr::all_of(time_var)) |>
    dplyr::mutate("{time_var}" := .sens_as_numeric_time(.data[[time_var]]))
  
  n_time_unique <- length(unique(df0[[time_var]]))
  deg_use <- min(poly_degree, max(1L, min(3L, n_time_unique - 1L)))
  
  rhs       <- .sens_poly_terms(time_var, deg_use)
  fixed_f   <- as.formula(sprintf("%s ~ %s", var_name, rhs))
  mixture_f <- as.formula(paste("~", rhs))
  
  set.seed(seed)
  
  minit_nr <- tryCatch(
    lcmm::hlme(
      fixed   = fixed_f,
      random  = ~ -1,
      idiag   = FALSE,
      nwg     = FALSE,
      subject = "id",
      ng      = 1,
      data    = df0
    ),
    error = function(e) e
  )
  
  if (inherits(minit_nr, "error")) {
    return(list(
      fit_ok = FALSE,
      model = NULL,
      pred_curve = NULL,
      error = paste("1-class initial model failed:", minit_nr$message)
    ))
  }
  
  model_K <- tryCatch(
    lcmm::gridsearch(
      rep     = 20,
      maxiter = 1000,
      minit   = minit_nr,
      hlme(
        fixed   = fixed_f,
        mixture = mixture_f,
        random  = ~ -1,
        idiag   = FALSE,
        nwg     = FALSE,
        subject = "id",
        ng      = K,
        data    = df0
      )
    ),
    error = function(e) e
  )
  
  if (inherits(model_K, "error")) {
    return(list(
      fit_ok = FALSE,
      model = NULL,
      pred_curve = NULL,
      error = paste("Gridsearch failed:", model_K$message)
    ))
  }
  
  pred_curve <- .sens_pred_curve_from_model(model_K, df0, time_var, K)
  
  if (is.null(pred_curve)) {
    return(list(
      fit_ok = FALSE,
      model = model_K,
      pred_curve = NULL,
      error = "No pred_curve could be extracted from fitted model."
    ))
  }
  
  list(
    fit_ok = TRUE,
    model = model_K,
    pred_curve = pred_curve,
    error = NULL
  )
}

# Compare a perturbed run's curve to the reference; returns tidy diagnostic
.sens_compare_to_reference <- function(fit_res, ref_curve, outcome,
                                       group_var, group_val, seed, K = 4) {
  
  if (isFALSE(fit_res$fit_ok) || is.null(fit_res$pred_curve)) {
    return(tibble::tibble(
      outcome    = outcome,
      !!group_var := group_val,
      seed       = seed,
      fit_ok     = FALSE,
      class_rank = NA_integer_,
      rmse       = NA_real_,
      maxabs     = NA_real_
    ))
  }
  
  aligned <- .sens_align_to_ref_rmse(fit_res$pred_curve, ref_curve, K = K)
  
  dplyr::inner_join(
    aligned, ref_curve,
    by = c("class_rank", "yearv"),
    suffix = c("", "_ref")
  ) |>
    dplyr::group_by(class_rank) |>
    dplyr::summarise(
      rmse   = sqrt(mean((pred - pred_ref)^2, na.rm = TRUE)),
      maxabs = max(abs(pred - pred_ref), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      outcome    = outcome,
      !!group_var := group_val,
      seed       = seed,
      fit_ok     = TRUE
    )
}

# Single-curve panel plot for visualisation
.sens_curve_panel_plot <- function(curve, outcome, title, K = 4) {
  
  ylab <- if (outcome == "dgcc") {
    "Effective glycaemic control coverage"
  } else {
    "Diabetic retinopathy screening coverage"
  }
  
  class_colors <- stats::setNames(ggsci::pal_lancet()(K), as.character(1:K))
  
  ggplot2::ggplot(
    curve,
    ggplot2::aes(x = yearv, y = pred,
                 colour = factor(class_rank, levels = 1:K))
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_colour_manual(
      values = class_colors,
      breaks = as.character(1:K),
      labels = paste("Class", 1:K),
      name   = NULL
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 12, by = 2),
      labels = seq(2011, 2023, by = 2)
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(title = title, x = "Year", y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 9),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(size = 8),
      axis.text        = ggplot2::element_text(size = 7),
      axis.title       = ggplot2::element_text(size = 10, face = "bold"),
      panel.grid       = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(colour = "black")
    )
}

# Assemble two-row figure (A = EGCC, B = DRSC) with separate boxed rows
.sens_assemble_two_row_figure <- function(panels_dgcc, panels_drsc) {
  
  strip_x_title <- ggplot2::theme(axis.title.x = ggplot2::element_blank())
  strip_y_all   <- ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y  = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )
  
  apply_strip <- function(pl) {
    purrr::imap(pl, function(p, i) {
      p <- p + strip_x_title
      if (i > 1) p <- p + strip_y_all
      p
    })
  }
  
  draw_box <- function() {
    cowplot::ggdraw() +
      cowplot::draw_line(x = c(0, 1), y = c(1, 1),
                         colour = "black", linewidth = 0.8) +
      cowplot::draw_line(x = c(0, 1), y = c(0, 0),
                         colour = "black", linewidth = 0.8) +
      cowplot::draw_line(x = c(0, 0), y = c(0, 1),
                         colour = "black", linewidth = 0.8) +
      cowplot::draw_line(x = c(1, 1), y = c(0, 1),
                         colour = "black", linewidth = 0.8)
  }
  
  grid_dgcc <- patchwork::wrap_plots(
    apply_strip(panels_dgcc), nrow = 1, guides = "collect"
  ) & ggplot2::theme(legend.position = "bottom")
  
  grid_drsc <- patchwork::wrap_plots(
    apply_strip(panels_drsc), nrow = 1, guides = "collect"
  ) & ggplot2::theme(legend.position = "bottom")
  
  row_A <- cowplot::ggdraw() +
    cowplot::draw_plot(grid_dgcc, x = 0.06, y = 0.02,
                       width = 0.92, height = 0.96) +
    cowplot::draw_plot(draw_box(), x = 0, y = 0, width = 1, height = 1) +
    cowplot::draw_label("A", x = 0.015, y = 0.985,
                        fontface = "bold", size = 16,
                        hjust = 0, vjust = 1)
  
  row_B <- cowplot::ggdraw() +
    cowplot::draw_plot(grid_drsc, x = 0.06, y = 0.02,
                       width = 0.92, height = 0.96) +
    cowplot::draw_plot(draw_box(), x = 0, y = 0, width = 1, height = 1) +
    cowplot::draw_label("B", x = 0.015, y = 0.985,
                        fontface = "bold", size = 16,
                        hjust = 0, vjust = 1)
  
  fig <- cowplot::plot_grid(row_A, row_B, ncol = 1,
                            rel_heights = c(1, 1))
  
  cowplot::ggdraw(fig) +
    cowplot::draw_label("Year", x = 0.52, y = 0.58, size = 10) +
    cowplot::draw_label("Year", x = 0.52, y = 0.08, size = 10)
}

# =========================================================================
# STAGE 1: Municipality subsampling validation
# =========================================================================
#
# Refits the selected 4-class cubic NRE model on random subsamples of
# municipalities (70%, 80%, 90% of IDs by default). Compares predicted
# class-mean trajectories to the full-sample reference via RMSE and
# max |Δ|, after class re-alignment.
#
# Returns: list with results (per seed × prop × class), summary table
# (median/IQR/max across seeds), and saved Figure 7.
# =========================================================================

make_subsampling_validation <- function(coverage_data,
                                        all_named_models,
                                        props     = c(0.90, 0.80, 0.70),
                                        seed_vec  = 1L,
                                        K         = 4L,
                                        time_var  = "year",
                                        out_path  = "outputs/figures/fig7_subsampling.png",
                                        width = 12, height = 10, dpi = 300) {
  
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  # ---- Build reference curves from saved full-sample models ----
  
  df_full <- coverage_data |>
    tidyr::drop_na(id, dgcc, drsc, dplyr::all_of(time_var)) |>
    dplyr::mutate("{time_var}" := .sens_as_numeric_time(.data[[time_var]]))
  
  ref_dgcc <- .sens_pred_curve_from_model(
    .sens_ref_model_for("dgcc", all_named_models),
    df_full, time_var = time_var, K = K
  )
  
  ref_drsc <- .sens_pred_curve_from_model(
    .sens_ref_model_for("drsc", all_named_models),
    df_full, time_var = time_var, K = K
  )
  
  # ---- Run subsamples across seed × prop × outcome ----
  
  results <- purrr::map_dfr(c("dgcc", "drsc"), function(var) {
    
    ref <- if (var == "dgcc") ref_dgcc else ref_drsc
    
    purrr::map_dfr(seed_vec, function(s) {
      
      purrr::map_dfr(props, function(p) {
        
        message(sprintf("[Stage 1] %s | prop = %.0f%% | seed = %d",
                        var, p * 100, s))
        
        set.seed(s)
        ids       <- unique(coverage_data$id)
        ids_sub   <- sample(ids, size = max(1L, floor(p * length(ids))))
        data_sub  <- dplyr::filter(coverage_data, id %in% ids_sub)
        
        fit <- .sens_fit_selected_model(
          var_name = var,
          data_in  = data_sub,
          K        = K,
          time_var = time_var,
          seed     = s
        )
        
        .sens_compare_to_reference(
          fit_res   = fit,
          ref_curve = ref,
          outcome   = var,
          group_var = "prop",
          group_val = p,
          seed      = s,
          K         = K
        )
      })
    })
  })
  
  # ---- Summary table across seeds ----
  
  summary_tbl <- results |>
    dplyr::filter(fit_ok) |>
    dplyr::group_by(outcome, prop, class_rank) |>
    dplyr::summarise(
      n_success     = dplyr::n(),
      median_rmse   = median(rmse,   na.rm = TRUE),
      iqr_rmse      = stats::IQR(rmse,   na.rm = TRUE),
      max_rmse      = max(rmse,      na.rm = TRUE),
      median_maxabs = median(maxabs, na.rm = TRUE),
      iqr_maxabs    = stats::IQR(maxabs, na.rm = TRUE),
      max_maxabs    = max(maxabs,    na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(outcome, prop, class_rank)
  
  # ---- Figure 7: representative seed panels ----
  
  plot_seed <- seed_vec[1]
  
  build_subsample_panels <- function(var, ref_curve) {
    
    orig_panel <- .sens_curve_panel_plot(ref_curve, var, "Original")
    
    prop_panels <- purrr::map(props, function(p) {
      
      set.seed(plot_seed)
      ids      <- unique(coverage_data$id)
      ids_sub  <- sample(ids, size = max(1L, floor(p * length(ids))))
      data_sub <- dplyr::filter(coverage_data, id %in% ids_sub)
      
      fit <- .sens_fit_selected_model(
        var_name = var,
        data_in  = data_sub,
        K        = K,
        time_var = time_var,
        seed     = plot_seed
      )
      
      curve <- if (!is.null(fit$pred_curve)) {
        .sens_align_to_ref_rmse(fit$pred_curve, ref_curve, K = K)
      } else {
        ref_curve
      }
      
      .sens_curve_panel_plot(curve, var,
                             paste0(p * 100, "% subsample"))
    })
    
    c(list(orig_panel), prop_panels)
  }
  
  panels_dgcc <- build_subsample_panels("dgcc", ref_dgcc)
  panels_drsc <- build_subsample_panels("drsc", ref_drsc)
  
  final_fig <- .sens_assemble_two_row_figure(panels_dgcc, panels_drsc)
  
  ggplot2::ggsave(out_path, plot = final_fig,
                  width = width, height = height,
                  dpi = dpi, units = "in")
  
  list(
    results   = results,
    summary   = summary_tbl,
    plot_path = out_path
  )
}


# =========================================================================
# STAGE 2: Municipality-year observation perturbation
# =========================================================================
#
# Refits the selected 4-class cubic NRE model on perturbed datasets
# where a fraction (5%, 10%, 15%, 20% by default) of municipality-year
# observations are randomly removed. Municipalities with fewer than
# `min_timepoints` remaining time points are excluded.
#
# Compares predicted class-mean trajectories to the full-sample
# reference via RMSE and max |Δ|, after class re-alignment.
#
# Returns: list with results (per seed × drop_prop × class), summary
# table, and saved Figure 8.
# =========================================================================

make_perturbation_validation <- function(coverage_data,
                                         all_named_models,
                                         drop_props     = c(0.05, 0.10, 0.15, 0.20),
                                         seed_vec       = 1L,
                                         K              = 4L,
                                         time_var       = "year",
                                         min_timepoints = 2L,
                                         out_path  = "outputs/figures/fig8_perturbation.png",
                                         width = 12, height = 10, dpi = 300) {
  
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  # ---- Build reference curves from saved full-sample models ----
  
  df_full <- coverage_data |>
    tidyr::drop_na(id, dgcc, drsc, dplyr::all_of(time_var)) |>
    dplyr::mutate("{time_var}" := .sens_as_numeric_time(.data[[time_var]]))
  
  ref_dgcc <- .sens_pred_curve_from_model(
    .sens_ref_model_for("dgcc", all_named_models),
    df_full, time_var = time_var, K = K
  )
  
  ref_drsc <- .sens_pred_curve_from_model(
    .sens_ref_model_for("drsc", all_named_models),
    df_full, time_var = time_var, K = K
  )
  
  # ---- Perturbation helper ----
  
  make_perturbed_data <- function(data, var, dp, seed) {
    
    set.seed(seed)
    
    df0 <- data |>
      tidyr::drop_na(id, dplyr::all_of(var), dplyr::all_of(time_var))
    
    n_drop   <- floor(dp * nrow(df0))
    drop_idx <- sample(seq_len(nrow(df0)), size = n_drop, replace = FALSE)
    
    df0[-drop_idx, ] |>
      dplyr::group_by(id) |>
      dplyr::filter(dplyr::n_distinct(.data[[time_var]]) >= min_timepoints) |>
      dplyr::ungroup()
  }
  
  # ---- Run perturbations across seed × drop_prop × outcome ----
  
  results <- purrr::map_dfr(c("dgcc", "drsc"), function(var) {
    
    ref <- if (var == "dgcc") ref_dgcc else ref_drsc
    
    purrr::map_dfr(seed_vec, function(s) {
      
      purrr::map_dfr(drop_props, function(dp) {
        
        message(sprintf("[Stage 2] %s | removed = %.0f%% | seed = %d",
                        var, dp * 100, s))
        
        data_pert <- make_perturbed_data(coverage_data, var, dp, s)
        
        fit <- .sens_fit_selected_model(
          var_name = var,
          data_in  = data_pert,
          K        = K,
          time_var = time_var,
          seed     = s
        )
        
        .sens_compare_to_reference(
          fit_res   = fit,
          ref_curve = ref,
          outcome   = var,
          group_var = "drop_prop",
          group_val = dp,
          seed      = s,
          K         = K
        )
      })
    })
  })
  
  # ---- Summary table across seeds ----
  
  summary_tbl <- results |>
    dplyr::filter(fit_ok) |>
    dplyr::group_by(outcome, drop_prop, class_rank) |>
    dplyr::summarise(
      n_success     = dplyr::n(),
      median_rmse   = median(rmse,   na.rm = TRUE),
      iqr_rmse      = stats::IQR(rmse,   na.rm = TRUE),
      max_rmse      = max(rmse,      na.rm = TRUE),
      median_maxabs = median(maxabs, na.rm = TRUE),
      iqr_maxabs    = stats::IQR(maxabs, na.rm = TRUE),
      max_maxabs    = max(maxabs,    na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(outcome, drop_prop, class_rank)
  
  # ---- Figure 8: representative seed panels ----
  
  plot_seed <- seed_vec[1]
  
  build_perturbation_panels <- function(var, ref_curve) {
    
    orig_panel <- .sens_curve_panel_plot(ref_curve, var, "Original")
    
    drop_panels <- purrr::map(drop_props, function(dp) {
      
      data_pert <- make_perturbed_data(coverage_data, var, dp, plot_seed)
      
      fit <- .sens_fit_selected_model(
        var_name = var,
        data_in  = data_pert,
        K        = K,
        time_var = time_var,
        seed     = plot_seed
      )
      
      curve <- if (!is.null(fit$pred_curve)) {
        .sens_align_to_ref_rmse(fit$pred_curve, ref_curve, K = K)
      } else {
        ref_curve
      }
      
      .sens_curve_panel_plot(curve, var,
                             paste0(dp * 100, "% removed"))
    })
    
    c(list(orig_panel), drop_panels)
  }
  
  panels_dgcc <- build_perturbation_panels("dgcc", ref_dgcc)
  panels_drsc <- build_perturbation_panels("drsc", ref_drsc)
  
  final_fig <- .sens_assemble_two_row_figure(panels_dgcc, panels_drsc)
  
  ggplot2::ggsave(out_path, plot = final_fig,
                  width = width, height = height,
                  dpi = dpi, units = "in")
  
  list(
    results   = results,
    summary   = summary_tbl,
    plot_path = out_path
  )
}
