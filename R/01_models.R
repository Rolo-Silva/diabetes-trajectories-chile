# R/01_models.R
# -------------------------------------------------------------------------
# Model adequacy, diagnostics, and visualisation for LCMM trajectory
# models. Reads pre-fitted models from all_named_models.rds (built via
# scripts/fit_models.R).
#
# Functions:
#   - make_model_adequacy_table():    full diagnostic table (BIC, entropy,
#                                     OCC, APPA, Lennon DoS, VLMR-LRT)
#   - make_bic_elbow_plot():          Figure 2 — BIC elbow per structure
#   - make_predicted_means_plot():    Figure 3 — predicted means A-J
#   - make_spaghetti_plot():          Figure 4 — observed + fitted (Struct D)
#   - make_elsensohn_plot():          Figure 5 — residual envelope (Struct D)
#   - make_residual_benchmark_plot(): Figure 6 — residuals (Struct A)
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# make_model_adequacy_table()
# -------------------------------------------------------------------------

make_model_adequacy_table <- function(all_named_models, coverage_data) {
  
  # --- Helper: build summary table from lcmm summarytable() ---------------
  
  build_summary_table <- function(model_names, model_list) {
    required_metrics <- c("G", "loglik", "conv", "npm", "AIC", "BIC",
                          "SABIC", "entropy", "ICL1", "ICL2", "%class")
    summary_list <- lapply(model_names, function(model_name) {
      model_obj <- model_list[[model_name]]
      if (is.null(model_obj)) return(tibble::tibble(Model = model_name))
      tryCatch({
        summ      <- lcmm::summarytable(model_obj, which = required_metrics)
        first_row <- as.data.frame(summ)[1, , drop = FALSE]
        for (i in 1:7) {
          if (!paste0("%class", i) %in% names(first_row))
            first_row[[paste0("%class", i)]] <- NA_real_
        }
        for (m in setdiff(required_metrics, names(first_row)))
          first_row[[m]] <- NA_real_
        first_row$Model <- model_name
        first_row[, c("Model", "G", "loglik", "conv", "npm", "AIC", "BIC",
                      "SABIC", "entropy", "ICL1", "ICL2", paste0("%class", 1:7))]
      }, error = function(e) tibble::tibble(Model = model_name))
    })
    dplyr::bind_rows(summary_list)
  }
  
  # --- Helper: smallest class size and count from postprob() --------------
  
  extract_postprob <- function(model) {
    tryCatch({
      pp <- lcmm::postprob(model)[[1]]
      list(size = min(pp[2, ], na.rm = TRUE), count = min(pp[1, ], na.rm = TRUE))
    }, error = function(e) list(size = NA, count = NA))
  }
  
  # --- Helper: OCC, APPA, mismatch from LCTMtoolkit ----------------------
  
  extract_occ_appa_mismatch <- function(model) {
    tryCatch({
      tk <- LCTMtools::LCTMtoolkit(model)
      list(
        occ      = min(as.numeric(tk$occ[1, ]),      na.rm = TRUE),
        appa     = min(as.numeric(tk$appa[1, ]),     na.rm = TRUE),
        mismatch = max(as.numeric(tk$mismatch[1, ]), na.rm = TRUE)
      )
    }, error = function(e) list(occ = NA, appa = NA, mismatch = NA))
  }
  
  # --- Helper: Vuong-Lo-Mendell-Rubin LRT p-value -------------------------
  
  extract_vllrt <- function(prev_model, curr_model) {
    tryCatch({
      if (is.null(prev_model) || is.null(curr_model)) return(NA_real_)
      lrt_stat <- 2 * (curr_model$loglik - prev_model$loglik)
      df <- lcmm::summarytable(curr_model)[1, "npm"] -
        lcmm::summarytable(prev_model)[1, "npm"]
      pchisq(lrt_stat, df = df, lower.tail = FALSE)
    }, error = function(e) NA_real_)
  }
  
  # --- Helper: combine all diagnostics per model --------------------------
  
  process_model <- function(model_name, model, prev_model = NULL) {
    tryCatch({
      if (is.null(model$ng)) return(tibble::tibble(Model = model_name, Error = "Invalid"))
      pp    <- if (model$ng > 1) extract_postprob(model) else list(size = NA, count = NA)
      occ   <- extract_occ_appa_mismatch(model)
      p_val <- if (!is.null(prev_model)) extract_vllrt(prev_model, model) else NA
      tibble::tibble(
        Model                          = model_name,
        Smallest_Class_Size_Percentage = pp$size,
        Smallest_Class_Count           = pp$count,
        Lowest_OCC                     = occ$occ,
        Lowest_APPA                    = occ$appa,
        Highest_Mismatch               = occ$mismatch,
        VLMRLRT_P_Value                = p_val
      )
    }, error = function(e) tibble::tibble(Model = model_name, Error = as.character(e)))
  }
  
  process_all_models <- function(model_list) {
    results <- list(); prev_model <- NULL
    for (i in seq_along(model_list)) {
      name        <- names(model_list)[i]
      results[[i]] <- process_model(name, model_list[[i]], prev_model)
      prev_model  <- model_list[[i]]
    }
    dplyr::bind_rows(results)
  }
  
  # --- Helper: Lennon-style DoS -------------------------------------------
  # Pairwise Mahalanobis distances between class mean trajectory vectors
  # (K x T matrix), weighted by class proportions, using the observed
  # outcome covariance matrix as the metric.
  
  calculate_lennon_dos <- function(model, data, outcome_var,
                                   time_var = "year", id_var = "id") {
    tryCatch({
      if (is.null(model) || is.null(model$ng)) return(NA_real_)
      K <- model$ng
      if (K <= 1) return(NA_real_)
      
      pred      <- model$pred
      pred_cols <- grep("^pred_m\\d+$", names(pred), value = TRUE)
      if (length(pred_cols) != K) return(NA_real_)
      
      data_a <- data |>
        tidyr::drop_na(dplyr::all_of(c(id_var, time_var, outcome_var))) |>
        dplyr::arrange(.data[[id_var]], .data[[time_var]]) |>
        dplyr::group_by(.data[[id_var]]) |>
        dplyr::mutate(.rid = dplyr::row_number()) |>
        dplyr::ungroup()
      
      pred_a <- pred |>
        dplyr::arrange(.data[[id_var]]) |>
        dplyr::group_by(.data[[id_var]]) |>
        dplyr::mutate(.rid = dplyr::row_number()) |>
        dplyr::ungroup()
      
      joined <- dplyr::left_join(
        dplyr::select(data_a, dplyr::all_of(c(id_var, time_var, outcome_var)), .rid),
        dplyr::select(pred_a, dplyr::all_of(id_var), .rid, dplyr::all_of(pred_cols)),
        by = c(id_var, ".rid")
      )
      
      M <- joined |>
        tidyr::pivot_longer(dplyr::all_of(pred_cols),
                            names_to = "cl", values_to = "pred") |>
        dplyr::mutate(cl = readr::parse_number(cl)) |>
        dplyr::group_by(cl, .data[[time_var]]) |>
        dplyr::summarise(pred = mean(pred, na.rm = TRUE), .groups = "drop") |>
        tidyr::pivot_wider(names_from = dplyr::all_of(time_var),
                           values_from = pred) |>
        dplyr::arrange(cl) |>
        dplyr::select(-cl) |>
        as.matrix()
      
      Y <- data_a |>
        dplyr::select(dplyr::all_of(c(id_var, time_var, outcome_var))) |>
        tidyr::pivot_wider(names_from = dplyr::all_of(time_var),
                           values_from = dplyr::all_of(outcome_var)) |>
        dplyr::select(-dplyr::all_of(id_var)) |>
        as.matrix()
      
      S     <- stats::cov(Y, use = "pairwise.complete.obs")
      S_inv <- tryCatch(solve(S), error = function(e) MASS::ginv(S))
      
      pp          <- suppressMessages(suppressWarnings(lcmm::postprob(model)[[1]]))
      class_props <- as.numeric(pp[2, ]) / 100
      
      dos <- 2 * sum(combn(seq_len(K), 2, function(ix) {
        d <- M[ix[1], ] - M[ix[2], ]
        sqrt(as.numeric(t(d) %*% S_inv %*% d)) *
          class_props[ix[1]] * class_props[ix[2]]
      }, simplify = TRUE), na.rm = TRUE)
      
      round(dos, 4)
    }, error = function(e) NA_real_)
  }
  
  add_lennon_dos <- function(model_list, table, coverage_data) {
    get_outcome <- function(name) {
      if (stringr::str_detect(name, "dgcc")) return("dgcc")
      if (stringr::str_detect(name, "drsc")) return("drsc")
      return(NA_character_)
    }
    
    df <- purrr::map_dfr(names(model_list), function(name) {
      outcome_var <- get_outcome(name)
      if (is.na(outcome_var)) return(tibble::tibble(Model = name, DoS = NA_real_))
      dos <- calculate_lennon_dos(model_list[[name]], coverage_data,
                                  outcome_var = outcome_var)
      tibble::tibble(Model = name, DoS = dos)
    })
    dplyr::left_join(table, df, by = "Model")
  }
  
  # --- Run the pipeline ---------------------------------------------------
  
  summary_table    <- build_summary_table(names(all_named_models), all_named_models)
  diagnostic_table <- process_all_models(all_named_models)
  
  summary_table$Model    <- as.character(summary_table$Model)
  diagnostic_table$Model <- as.character(diagnostic_table$Model)
  
  adequacy_table <- dplyr::left_join(summary_table, diagnostic_table, by = "Model")
  final_table    <- add_lennon_dos(all_named_models, adequacy_table, coverage_data)
  
  class_cols <- grep("^%class\\d+$", names(final_table), value = TRUE)
  
  model_adequacy_table <- final_table |>
    dplyr::mutate(structure = dplyr::case_when(
      stringr::str_detect(Model, "linear_nre_homocedastic")       ~ "A",
      stringr::str_detect(Model, "linear_nre_heterocedastic")     ~ "B",
      stringr::str_detect(Model, "quadratic_nre")                 ~ "C",
      stringr::str_detect(Model, "cubic_nre")                     ~ "D",
      stringr::str_detect(Model, "linear_random_intercept_slope") ~ "F",
      stringr::str_detect(Model, "linear_random_intercept")       ~ "E",
      stringr::str_detect(Model, "quadratic_random_effects_prop") ~ "H",
      stringr::str_detect(Model, "quadratic_random_effects")      ~ "G",
      stringr::str_detect(Model, "cubic_random_effects_prop")     ~ "J",
      stringr::str_detect(Model, "cubic_random_effects")          ~ "I",
      TRUE ~ "Other"
    )) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      Class_Sizes_Percentage = if (is.na(G) || G <= 1) NA_character_ else {
        vals <- dplyr::c_across(dplyr::all_of(class_cols))
        vals <- vals[!is.na(vals)]
        paste0(sprintf("%.1f", vals), collapse = "; ")
      }
    ) |>
    dplyr::ungroup()
  
  return(model_adequacy_table)
}


# -------------------------------------------------------------------------
# make_bic_elbow_plot()
# -------------------------------------------------------------------------

make_bic_elbow_plot <- function(model_adequacy_table,
                                out_path = "outputs/figures/fig2_bic_elbow.png",
                                width = 12, height = 8, dpi = 300) {
  
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  elbow_tbl <- model_adequacy_table
  names(elbow_tbl) <- tolower(names(elbow_tbl))
  
  if ("model_name" %in% names(elbow_tbl)) elbow_tbl <- dplyr::rename(elbow_tbl, model = model_name)
  if ("modelo"     %in% names(elbow_tbl)) elbow_tbl <- dplyr::rename(elbow_tbl, model = modelo)
  if ("bic_value"  %in% names(elbow_tbl)) elbow_tbl <- dplyr::rename(elbow_tbl, bic   = bic_value)
  stopifnot(all(c("model", "bic") %in% names(elbow_tbl)))
  
  elbow_tbl <- elbow_tbl |>
    dplyr::mutate(
      Model = as.character(model),
      BIC   = if (is.numeric(bic)) bic else readr::parse_number(as.character(bic))
    ) |>
    dplyr::select(-dplyr::any_of("bic_display")) |>
    dplyr::filter(is.finite(BIC))
  
  g_col <- names(elbow_tbl)[tolower(names(elbow_tbl)) %in% c("g", "k", "classes", "n_classes", "nclass")]
  elbow_tbl <- elbow_tbl |>
    dplyr::mutate(
      G = if (length(g_col)) suppressWarnings(as.integer(.data[[g_col[1]]]))
      else suppressWarnings(as.integer(stringr::str_extract(Model, "^[0-9]+")))
    )
  
  outcome_col <- names(elbow_tbl)[tolower(names(elbow_tbl)) %in% c("outcome", "variable", "outcome_name")]
  elbow_tbl <- elbow_tbl |>
    dplyr::mutate(
      Outcome = if (length(outcome_col)) toupper(as.character(.data[[outcome_col[1]]]))
      else {
        out <- stringr::str_extract(Model, "(dgcc|drsc)")
        toupper(ifelse(is.na(out), "UNKNOWN", out))
      }
    )
  
  elbow_tbl <- elbow_tbl |>
    dplyr::mutate(
      Structure = Model |>
        stringr::str_replace("^[0-9]+class_", "") |>
        stringr::str_replace("_model.*$", "") |>
        stringr::str_replace("_(dgcc|drsc)$", ""),
      Outcome = ifelse(Outcome == "DGCC", "EGCC", Outcome)
    ) |>
    dplyr::filter(!is.na(G), G > 0, Outcome %in% c("EGCC", "DRSC"), !is.na(Structure))
  
  structure_labels <- elbow_tbl |>
    dplyr::mutate(
      Label = dplyr::case_when(
        stringr::str_detect(Structure, "linear_nre_homocedastic")       ~ "A",
        stringr::str_detect(Structure, "linear_nre_heterocedastic")     ~ "B",
        stringr::str_detect(Structure, "quadratic_nre")                 ~ "C",
        stringr::str_detect(Structure, "cubic_nre")                     ~ "D",
        stringr::str_detect(Structure, "linear_random_intercept$")      ~ "E",
        stringr::str_detect(Structure, "linear_random_intercept_slope") ~ "F",
        stringr::str_detect(Structure, "quadratic_random_effects$")     ~ "G",
        stringr::str_detect(Structure, "quadratic_random_effects_prop") ~ "H",
        stringr::str_detect(Structure, "cubic_random_effects$")         ~ "I",
        stringr::str_detect(Structure, "cubic_random_effects_prop")     ~ "J",
        TRUE ~ "Other"
      ),
      PrettyTitle = paste("Structure", Label, "|", Outcome)
    ) |>
    dplyr::distinct(Structure, Outcome, PrettyTitle)
  
  prop_cols <- names(elbow_tbl)[stringr::str_detect(names(elbow_tbl), "^%class")]
  if (!length(prop_cols))
    prop_cols <- names(elbow_tbl)[stringr::str_detect(names(elbow_tbl), "^class\\d+")]
  stopifnot(length(prop_cols) > 0)
  
  class_df_long <- elbow_tbl |>
    tidyr::pivot_longer(dplyr::all_of(prop_cols), names_to = "Class", values_to = "Proportion") |>
    dplyr::mutate(
      Class      = factor(Class, levels = prop_cols),
      Proportion = suppressWarnings(as.numeric(Proportion))
    )
  
  valid_keys <- class_df_long |>
    dplyr::distinct(Outcome, Structure, G) |>
    dplyr::count(Outcome, Structure, name = "nG") |>
    dplyr::filter(nG >= 2L) |>
    dplyr::select(Outcome, Structure)
  
  plot_keys <- valid_keys |>
    dplyr::left_join(structure_labels, by = c("Outcome", "Structure"))
  
  plot_elbow_stack <- function(df_model, title = NULL) {
    df_model <- df_model |> dplyr::filter(is.finite(BIC))
    
    df_summary <- df_model |>
      dplyr::group_by(G, Class) |>
      dplyr::summarise(Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop")
    
    bic_df <- df_model |>
      dplyr::group_by(G) |>
      dplyr::summarise(BIC = mean(BIC, na.rm = TRUE), .groups = "drop")
    
    if (nrow(bic_df) < 2) return(NULL)
    
    bmin <- min(bic_df$BIC, na.rm = TRUE)
    bmax <- max(bic_df$BIC, na.rm = TRUE)
    sf   <- ifelse(bmax > bmin, 100 / (bmax - bmin), 1)
    bic_df <- bic_df |> dplyr::mutate(BIC_scaled = (BIC - bmin) * sf)
    
    ggplot2::ggplot(df_summary, ggplot2::aes(x = factor(G), y = Proportion, fill = Class)) +
      ggplot2::geom_col(colour = "white", linewidth = 0.15) +
      ggplot2::scale_fill_brewer(palette = "Reds") +
      ggplot2::geom_line(
        data = bic_df,
        ggplot2::aes(x = factor(G), y = BIC_scaled, group = 1),
        inherit.aes = FALSE, colour = "black", linewidth = 0.45
      ) +
      ggplot2::geom_point(
        data = bic_df,
        ggplot2::aes(x = factor(G), y = BIC_scaled),
        inherit.aes = FALSE, colour = "black", size = 0.9
      ) +
      ggplot2::scale_y_continuous(
        name      = "% of municipalities",
        limits    = c(0, 100),
        sec.axis  = ggplot2::sec_axis(~ . / sf + bmin, name = "BIC")
      ) +
      ggplot2::labs(title = title, x = "Number of classes", fill = "Class") +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        legend.position  = "none",
        plot.title       = ggplot2::element_text(size = 7.5, face = "bold", hjust = 0.5),
        axis.title       = ggplot2::element_blank()
      )
  }
  
  plot_data_list <- purrr::map2(
    plot_keys$Outcome, plot_keys$Structure,
    ~ class_df_long |> dplyr::filter(Outcome == .x, Structure == .y)
  )
  plots <- purrr::map2(plot_data_list, plot_keys$PrettyTitle, plot_elbow_stack)
  valid_idx <- purrr::map_lgl(plots, ~ !is.null(.))
  plots     <- plots[valid_idx]
  names(plots) <- plot_keys$PrettyTitle[valid_idx]
  
  extract_label <- function(name) stringr::str_remove(stringr::str_extract(name, "Structure [A-J]"), "Structure ")
  
  order_by_label <- function(pl) {
    if (!length(pl)) return(pl)
    ord <- tibble::tibble(name = names(pl)) |>
      dplyr::mutate(label = extract_label(name),
                    order = match(label, LETTERS[1:10])) |>
      dplyr::arrange(order)
    pl[ord$name]
  }
  
  remove_inner_axes <- function(pl, ncol = 5) {
    n          <- length(pl)
    nrow       <- ceiling(n / ncol)
    left_idx   <- seq(1, by = ncol, length.out = nrow)
    bottom_idx <- ((nrow - 1) * ncol + 1):min(nrow * ncol, n)
    right_idx  <- pmin(seq(ncol, by = ncol, length.out = nrow), n)
    
    purrr::imap(pl, function(p, i) {
      p <- p + ggplot2::theme(
        axis.text.x       = ggplot2::element_blank(),
        axis.ticks.x      = ggplot2::element_blank(),
        axis.text.y       = ggplot2::element_blank(),
        axis.ticks.y      = ggplot2::element_blank(),
        axis.text.y.right = ggplot2::element_blank(),
        axis.ticks.y.right = ggplot2::element_blank()
      )
      if (i %in% left_idx)
        p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(),
                                axis.ticks.y = ggplot2::element_line())
      if (i %in% bottom_idx)
        p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(),
                                axis.ticks.x = ggplot2::element_line())
      if (i %in% right_idx)
        p <- p + ggplot2::theme(axis.text.y.right = ggplot2::element_text(),
                                axis.ticks.y.right = ggplot2::element_line())
      p
    })
  }
  
  add_panel_labels <- function(p, tag) {
    cowplot::ggdraw(p) +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(colour = "grey55", fill = NA, linewidth = 0.5),
        plot.margin     = ggplot2::margin(24, 34, 30, 34)
      ) +
      cowplot::draw_label(tag,                  x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 13, fontface = "bold") +
      cowplot::draw_label("% of municipalities", x = 0.01, y = 0.50, angle = 90, size = 9) +
      cowplot::draw_label("Number of classes",   x = 0.50, y = 0.02, size = 9) +
      cowplot::draw_label("BIC",                 x = 0.99, y = 0.50, angle = 90, size = 9)
  }
  
  egcc_list <- plots[grepl("\\|\\s*EGCC$", names(plots))]
  drsc_list <- plots[grepl("\\|\\s*DRSC$", names(plots))]
  
  egcc_list <- utils::head(order_by_label(egcc_list), 10)
  drsc_list <- utils::head(order_by_label(drsc_list), 10)
  
  egcc_grid <- patchwork::wrap_plots(remove_inner_axes(egcc_list, ncol = 5), ncol = 5)
  drsc_grid <- patchwork::wrap_plots(remove_inner_axes(drsc_list, ncol = 5), ncol = 5)
  
  p_top    <- add_panel_labels(egcc_grid, "A")
  p_bottom <- add_panel_labels(drsc_grid, "B")
  
  final_fig <- p_top / p_bottom
  
  ggplot2::ggsave(out_path, plot = final_fig, width = width, height = height,
                  dpi = dpi, units = "in")
  
  return(out_path)
}


# -------------------------------------------------------------------------
# make_predicted_means_plot()
# -------------------------------------------------------------------------

make_predicted_means_plot <- function(all_named_models,
                                      out_path = "outputs/figures/fig3_predicted_means.png",
                                      width = 12, height = 8, dpi = 300) {
  
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  structure_label <- function(model_name) {
    dplyr::case_when(
      stringr::str_detect(model_name, "linear_nre_homocedastic")       ~ "A",
      stringr::str_detect(model_name, "linear_nre_heterocedastic")     ~ "B",
      stringr::str_detect(model_name, "quadratic_nre")                 ~ "C",
      stringr::str_detect(model_name, "cubic_nre")                     ~ "D",
      stringr::str_detect(model_name, "linear_random_intercept$")      ~ "E",
      stringr::str_detect(model_name, "linear_random_intercept_slope") ~ "F",
      stringr::str_detect(model_name, "quadratic_random_effects$")     ~ "G",
      stringr::str_detect(model_name, "quadratic_random_effects_prop") ~ "H",
      stringr::str_detect(model_name, "cubic_random_effects$")         ~ "I",
      stringr::str_detect(model_name, "cubic_random_effects_prop")     ~ "J",
      TRUE ~ "Other"
    )
  }
  
  model_names_4class <- names(all_named_models)[
    stringr::str_detect(names(all_named_models), "^4class_")
  ]
  
  build_panel <- function(model_name) {
    model <- all_named_models[[model_name]]
    outcome_raw <- stringr::str_extract(model_name, "dgcc|drsc")
    structure    <- structure_label(model_name)
    panel_title  <- paste0("Structure ", structure)
    
    pred      <- model$pred
    pred_cols <- grep("^pred_m[1-4]$", names(pred), value = TRUE)
    
    if (length(pred_cols) == 0) {
      return(
        ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
          ggplot2::geom_text(label = paste0("Structure ", structure, "\ndid not converge"),
                             size = 3, hjust = 0.5, vjust = 0.5) +
          ggplot2::xlim(0, 2) + ggplot2::ylim(0, 2) +
          ggplot2::theme_void() +
          ggplot2::ggtitle(panel_title) +
          ggplot2::theme(
            plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold", size = 7.5),
            panel.border = ggplot2::element_rect(color = "grey80", fill = NA)
          )
      )
    }
    
    data <- eval(model$call$data)
    
    data <- data |>
      dplyr::group_by(id) |>
      dplyr::mutate(row_id = dplyr::row_number()) |>
      dplyr::ungroup()
    
    pred <- pred |>
      dplyr::group_by(id) |>
      dplyr::mutate(row_id = dplyr::row_number()) |>
      dplyr::ungroup()
    
    combined <- dplyr::left_join(
      dplyr::select(data, id, row_id, year),
      dplyr::select(pred, id, row_id, dplyr::all_of(pred_cols)),
      by = c("id", "row_id")
    )
    
    pred_long <- combined |>
      tidyr::pivot_longer(dplyr::all_of(pred_cols),
                          names_to  = "class_label",
                          values_to = "pred") |>
      dplyr::mutate(
        class = as.integer(gsub("pred_m", "", class_label)),
        label = paste("Class", class)
      )
    
    pred_summary <- pred_long |>
      dplyr::group_by(class, label, year) |>
      dplyr::summarise(pred = mean(pred, na.rm = TRUE), .groups = "drop")
    
    ggplot2::ggplot(pred_summary,
                    ggplot2::aes(x = year, y = pred,
                                 color = label, group = label)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggsci::scale_color_lancet(name = NULL) +
      ggplot2::scale_x_continuous(breaks = seq(0, 12, 4),
                                  labels = seq(2011, 2023, 4)) +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::ggtitle(panel_title) +
      ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
      ggplot2::theme_bw(base_size = 7) +
      ggplot2::theme(
        legend.position   = "bottom",
        legend.direction  = "horizontal",
        legend.key.height = ggplot2::unit(6, "pt"),
        legend.key.width  = ggplot2::unit(10, "pt"),
        legend.text       = ggplot2::element_text(size = 6),
        legend.margin     = ggplot2::margin(t = 2, b = -2),
        legend.spacing.y  = ggplot2::unit(0, "pt"),
        axis.text         = ggplot2::element_text(size = 6),
        panel.grid        = ggplot2::element_blank(),
        panel.border      = ggplot2::element_blank(),
        axis.line         = ggplot2::element_line(color = "black"),
        axis.title        = ggplot2::element_blank(),
        plot.title        = ggplot2::element_text(hjust = 0.5, face = "bold", size = 7.5)
      )
  }
  
  all_panels <- purrr::set_names(
    purrr::map(model_names_4class, build_panel),
    model_names_4class
  )
  
  order_panels <- function(panels) {
    labels <- structure_label(names(panels))
    ord    <- match(labels, LETTERS[1:10])
    panels[order(ord)]
  }
  
  dgcc_panels <- order_panels(all_panels[stringr::str_detect(names(all_panels), "dgcc")])
  drsc_panels <- order_panels(all_panels[stringr::str_detect(names(all_panels), "drsc")])
  
  while (length(dgcc_panels) < 10) dgcc_panels[[length(dgcc_panels) + 1]] <- patchwork::plot_spacer()
  while (length(drsc_panels) < 10) drsc_panels[[length(drsc_panels) + 1]] <- patchwork::plot_spacer()
  
  names(dgcc_panels) <- paste0("Structure ", LETTERS[1:10])
  names(drsc_panels) <- paste0("Structure ", LETTERS[1:10])
  
  grid_A <- patchwork::wrap_plots(dgcc_panels, ncol = 5, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  
  grid_B <- patchwork::wrap_plots(drsc_panels, ncol = 5, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  
  add_panel_labels <- function(p, tag, left_lab) {
    cowplot::ggdraw(p) +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(colour = "grey60", fill = NA, linewidth = 0.8),
        plot.margin     = ggplot2::margin(24, 36, 18, 44)
      ) +
      cowplot::draw_label(tag,      x = 0.012, y = 0.986, hjust = 0, vjust = 1,
                          fontface = "bold", size = 14) +
      cowplot::draw_label(left_lab, x = 0,     y = 0.5,   angle = 90,
                          hjust = 0.5, vjust = 0.5, size = 10) +
      cowplot::draw_label("Year",   x = 0.5,   y = 0.115,
                          vjust = 0, hjust = 0.5, size = 10)
  }
  
  panel_A <- add_panel_labels(grid_A, "A", "Effective glycaemic control coverage")
  panel_B <- add_panel_labels(grid_B, "B", "Diabetic retinopathy screening coverage")
  
  final_fig <- panel_A / panel_B
  
  ggplot2::ggsave(out_path, plot = final_fig, width = width, height = height,
                  dpi = dpi, units = "in")
  
  return(out_path)
}


# -------------------------------------------------------------------------
# make_spaghetti_plot()
# -------------------------------------------------------------------------

make_spaghetti_plot <- function(all_named_models,
                                trajectory_labels,
                                out_path = "outputs/figures/fig4_spaghetti.png",
                                width = 12, height = 7, dpi = 300) {
  
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  model_names <- c("4class_cubic_nre_dgcc_model", "4class_cubic_nre_drsc_model")
  
  build_panel <- function(model_name) {
    model <- all_named_models[[model_name]]
    
    if (is.null(model) || is.null(model$pred) || is.null(model$pprob)) {
      message("Skipping: ", model_name)
      return(NULL)
    }
    
    outcome_var <- ifelse(stringr::str_detect(model_name, "dgcc"), "dgcc", "drsc")
    y_label     <- ifelse(outcome_var == "dgcc",
                          "Effective glycaemic control coverage",
                          "Diabetic retinopathy screening coverage")
    
    model_data <- eval(model$call$data) |>
      dplyr::arrange(id, year) |>
      dplyr::group_by(id) |>
      dplyr::mutate(row_id = dplyr::row_number()) |>
      dplyr::ungroup()
    
    pred_data <- model$pred |>
      dplyr::arrange(id) |>
      dplyr::group_by(id) |>
      dplyr::mutate(row_id = dplyr::row_number()) |>
      dplyr::ungroup()
    
    merged_df <- dplyr::left_join(model_data, pred_data, by = c("id", "row_id")) |>
      dplyr::left_join(model$pprob, by = "id") |>
      dplyr::filter(!is.na(class)) |>
      dplyr::mutate(class = as.integer(as.character(class)))
    
    class_counts <- merged_df |>
      dplyr::distinct(id, class) |>
      dplyr::group_by(class) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
      dplyr::mutate(percentage = n / sum(n))
    
    labels_vec <- trajectory_labels[[model_name]]
    label_df <- data.frame(
      class      = 1:4,
      traj_label = labels_vec,
      stringsAsFactors = FALSE
    ) |>
      dplyr::left_join(class_counts, by = "class") |>
      dplyr::mutate(
        label = paste0(traj_label,
                       " (", scales::percent(percentage, accuracy = 0.1), ")")
      )
    
    merged_df <- merged_df |>
      dplyr::left_join(dplyr::select(label_df, class, label), by = "class") |>
      dplyr::mutate(label = factor(label, levels = label_df$label))
    
    ggplot2::ggplot(merged_df,
                    ggplot2::aes(x = year, y = .data[[outcome_var]],
                                 group = id, color = label)) +
      ggplot2::geom_line(alpha = 0.15, linewidth = 0.3) +
      ggplot2::stat_summary(fun = mean, geom = "line",
                            ggplot2::aes(group = label),
                            linewidth = 0.9) +
      ggplot2::facet_wrap(~ label, nrow = 1, strip.position = "top") +
      ggsci::scale_color_lancet() +
      ggplot2::scale_y_continuous(limits = c(0, 1),
                                  labels = scales::percent) +
      ggplot2::scale_x_continuous(breaks = seq(0, 12, by = 4),
                                  labels = seq(2011, 2023, by = 4)) +
      ggplot2::labs(x = "Year", y = y_label) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        strip.text       = ggplot2::element_text(size = 8, face = "bold"),
        axis.text        = ggplot2::element_text(size = 7),
        axis.title       = ggplot2::element_text(size = 8.5),
        legend.position  = "none",
        panel.grid       = ggplot2::element_blank(),
        panel.border     = ggplot2::element_blank(),
        axis.line        = ggplot2::element_line(color = "black")
      )
  }
  
  panels <- purrr::compact(purrr::map(model_names, build_panel))
  stopifnot(length(panels) == 2)
  
  final_fig <- (panels[[1]] / panels[[2]]) +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(
      plot.tag          = ggplot2::element_text(size = 14, face = "bold"),
      plot.tag.position = c(0.015, 0.985),
      plot.background   = ggplot2::element_rect(colour = "grey70",
                                                fill = NA, linewidth = 0.6),
      plot.margin       = ggplot2::margin(10, 12, 12, 12)
    )
  
  ggplot2::ggsave(out_path, plot = final_fig, width = width, height = height,
                  dpi = dpi, units = "in")
  
  return(out_path)
}


# -------------------------------------------------------------------------
# make_elsensohn_plot()
# -------------------------------------------------------------------------

make_elsensohn_plot <- function(all_named_models,
                                trajectory_labels,
                                out_path = "outputs/figures/fig5_elsensohn.png",
                                width = 12, height = 5, dpi = 300) {
  
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  build_panel <- function(model_name, outcome_var, y_label) {
    model <- all_named_models[[model_name]]
    
    pred <- model$pred |>
      dplyr::group_by(id) |>
      dplyr::mutate(year = dplyr::row_number() - 1) |>
      dplyr::ungroup() |>
      dplyr::left_join(
        model$pprob |> dplyr::select(id, dplyr::starts_with("prob")),
        by = "id"
      )
    
    env_long <- pred |>
      dplyr::select(id, year, obs,
                    pred_m1, pred_m2, pred_m3, pred_m4,
                    prob1, prob2, prob3, prob4) |>
      tidyr::pivot_longer(
        cols = c(pred_m1, pred_m2, pred_m3, pred_m4),
        names_to = "pred_class", values_to = "mu_hat"
      ) |>
      dplyr::mutate(
        class = as.integer(gsub("pred_m", "", pred_class)),
        post_prob = dplyr::case_when(
          class == 1 ~ prob1, class == 2 ~ prob2,
          class == 3 ~ prob3, class == 4 ~ prob4
        ),
        resid = obs - mu_hat
      )
    
    env_summary <- env_long |>
      dplyr::group_by(year, class) |>
      dplyr::summarise(
        mu_hat    = dplyr::first(mu_hat),
        w_sum     = sum(post_prob, na.rm = TRUE),
        local_var = sum(post_prob * resid^2, na.rm = TRUE) / w_sum,
        local_sd  = sqrt(local_var),
        .groups   = "drop"
      )
    
    obs_traj <- pred |>
      dplyr::left_join(
        model$pprob |> dplyr::select(id, class),
        by = "id"
      ) |>
      dplyr::group_by(year, class) |>
      dplyr::summarise(obs_mean = mean(obs, na.rm = TRUE), .groups = "drop")
    
    env_plot <- env_summary |>
      dplyr::left_join(obs_traj, by = c("year", "class")) |>
      dplyr::mutate(
        upper_env_obs = obs_mean + local_sd,
        lower_env_obs = obs_mean - local_sd
      )
    
    class_sizes <- model$pprob |>
      dplyr::count(class) |>
      dplyr::mutate(pct = round(100 * n / sum(n), 1))
    
    labels_vec <- trajectory_labels[[model_name]]
    label_df <- data.frame(
      class      = 1:4,
      label_base = labels_vec,
      stringsAsFactors = FALSE
    ) |>
      dplyr::left_join(class_sizes, by = "class") |>
      dplyr::mutate(
        label = paste0(label_base, ": n=", n, " (", pct, "%)")
      )
    
    env_plot <- env_plot |>
      dplyr::left_join(dplyr::select(label_df, class, label), by = "class") |>
      dplyr::mutate(class_f = factor(label, levels = label_df$label))
    
    ggplot2::ggplot(
      env_plot,
      ggplot2::aes(x = year, y = obs_mean, colour = class_f, group = class_f)
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower_env_obs, ymax = upper_env_obs, fill = class_f),
        alpha = 0.10, colour = NA
      ) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::geom_line(ggplot2::aes(y = upper_env_obs),
                         linetype = "dashed", linewidth = 0.7) +
      ggplot2::geom_line(ggplot2::aes(y = lower_env_obs),
                         linetype = "dashed", linewidth = 0.7) +
      ggplot2::scale_x_continuous(breaks = 0:12, labels = 2011:2023) +
      ggplot2::scale_y_continuous(limits = c(0, 1),
                                  breaks = seq(0, 1, 0.2),
                                  labels = scales::percent) +
      ggsci::scale_color_lancet() +
      ggsci::scale_fill_lancet(guide = "none") +
      ggplot2::labs(x = "Year", y = y_label,
                    colour = "Latent Class", fill = "Latent Class") +
      ggplot2::guides(
        colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE,
                                       title.position = "top", title.hjust = 0),
        fill = "none"
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        legend.position      = c(0.5, 0.98),
        legend.justification = c(0.5, 1),
        legend.direction     = "horizontal",
        legend.background    = ggplot2::element_rect(
          fill = scales::alpha("white", 0.7), color = NA),
        legend.box.background = ggplot2::element_blank(),
        legend.margin        = ggplot2::margin(0, 0, 0, 0),
        legend.key.height    = grid::unit(8, "pt"),
        legend.key.width     = grid::unit(16, "pt"),
        panel.grid           = ggplot2::element_blank(),
        panel.border         = ggplot2::element_blank(),
        axis.line            = ggplot2::element_line(color = "black")
      )
  }
  
  panel_A <- build_panel(
    "4class_cubic_nre_dgcc_model", "dgcc",
    "Effective glycaemic control coverage"
  ) + ggplot2::labs(tag = "A")
  
  panel_B <- build_panel(
    "4class_cubic_nre_drsc_model", "drsc",
    "Diabetic retinopathy screening coverage"
  ) + ggplot2::labs(tag = "B")
  
  final_fig <- (panel_A | panel_B) &
    ggplot2::theme(
      plot.tag          = ggplot2::element_text(size = 18, face = "bold"),
      plot.tag.position = c(0.01, 0.98)
    )
  
  ggplot2::ggsave(out_path, plot = final_fig, width = width, height = height,
                  dpi = dpi, units = "in")
  
  return(out_path)
}


# -------------------------------------------------------------------------
# make_residual_benchmark_plot()
# -------------------------------------------------------------------------

make_residual_benchmark_plot <- function(all_named_models,
                                         out_path = "outputs/figures/fig6_residual_benchmark.png",
                                         width = 12, height = 7, dpi = 300) {
  
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  extract_residual_data <- function(model, outcome_var, model_label) {
    k     <- ifelse(is.null(model$ng), 1, model$ng)
    preds <- model$pred
    
    if (k == 1 && !"pred_ss1" %in% names(preds)) {
      if ("pred" %in% names(preds)) {
        preds$pred_ss1 <- preds$pred
      } else if ("pred_marg" %in% names(preds)) {
        preds$pred_ss1 <- preds$pred_marg
      } else {
        stop("Cannot find prediction column in model$pred.")
      }
    }
    
    nameofid   <- names(model$pred)[1]
    model_data <- eval(model$call$data)
    
    time_var <- intersect(c("ano", "year", "time", "Time"), names(model_data))
    if (length(time_var) == 0) stop("Cannot find time variable in model data.")
    time_var <- time_var[1]
    
    model_data <- model_data |>
      dplyr::group_by(.data[[nameofid]]) |>
      dplyr::mutate(row_id = dplyr::row_number()) |>
      dplyr::ungroup()
    
    preds <- preds |>
      dplyr::group_by(.data[[nameofid]]) |>
      dplyr::mutate(row_id = dplyr::row_number()) |>
      dplyr::ungroup()
    
    test <- dplyr::left_join(preds, model$pprob, by = nameofid) |>
      dplyr::left_join(
        dplyr::select(model_data, dplyr::all_of(nameofid),
                      dplyr::all_of(outcome_var),
                      dplyr::all_of(time_var),
                      row_id),
        by = c(nameofid, "row_id")
      )
    
    purrr::map_dfr(seq_len(k), function(i) {
      class_col <- if (k == 1) "pred_ss1" else paste0("pred_ss", i)
      test |>
        dplyr::filter(k == 1 | class == i) |>
        dplyr::transmute(
          ano_real    = .data[[time_var]],
          Residuals   = .data[[outcome_var]] - .data[[class_col]],
          class_label = paste("Class", i),
          model_label = model_label
        )
    })
  }
  
  residual_data <- dplyr::bind_rows(
    extract_residual_data(
      all_named_models[["4class_linear_nre_homocedastic_dgcc_model"]],
      "dgcc", "EGCC"
    ),
    extract_residual_data(
      all_named_models[["4class_linear_nre_homocedastic_drsc_model"]],
      "drsc", "DRSC"
    )
  )
  
  plot_residual_block <- function(df, title) {
    df <- df |> dplyr::filter(!is.na(class_label))
    
    ord_nums <- readr::parse_number(as.character(df$class_label))
    ord      <- paste("Class", sort(unique(ord_nums)))
    df       <- df |>
      dplyr::mutate(class_label = factor(class_label, levels = ord, ordered = TRUE))
    n_cls <- length(levels(df$class_label))
    
    ggplot2::ggplot(df, ggplot2::aes(x = ano_real, y = Residuals)) +
      ggplot2::geom_point(size = 0.2, alpha = 0.3) +
      ggplot2::stat_summary(fun = mean, geom = "line",
                            color = "darkcyan", linewidth = 0.6) +
      ggplot2::facet_wrap(~ class_label, nrow = 1, ncol = n_cls) +
      ggplot2::ylim(-1, 1) +
      ggplot2::scale_x_continuous(breaks = seq(2011, 2023, 2)) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(x = "Year", y = "Residuals", title = title) +
      ggplot2::theme(
        strip.text       = ggplot2::element_text(face = "bold", size = 8),
        axis.text.x      = ggplot2::element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y      = ggplot2::element_text(size = 7),
        plot.title       = ggplot2::element_text(face = "bold", size = 10, hjust = 0),
        panel.grid       = ggplot2::element_blank()
      )
  }
  
  panel_A <- plot_residual_block(
    residual_data |> dplyr::filter(model_label == "EGCC"), "A  EGCC"
  )
  panel_B <- plot_residual_block(
    residual_data |> dplyr::filter(model_label == "DRSC"), "B  DRSC"
  )
  
  final_fig <- (panel_A / panel_B)
  
  ggplot2::ggsave(out_path, plot = final_fig, width = width, height = height,
                  dpi = dpi, units = "in")
  
  return(out_path)
}