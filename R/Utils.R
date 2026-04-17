# R/Utils.R
# All pipeline functions for the trajectories targets project.
# Model fitting is handled separately in scripts/fit_models.R


# -------------------------------------------------------------------------
# make_all_series()
# Combines a list of yearly data frames and applies final cleaning/type casting.
# -------------------------------------------------------------------------

make_all_series <- function(series_list) {

  # Combine — safe because all inputs are character columns
  all_series <- purrr::list_rbind(series_list)

  # Initial type inference across the full combined data
  all_series <- type.convert(all_series, as.is = TRUE)

  # Standardise id_comuna to zero-padded 5-digit string
  all_series <- all_series |>
    dplyr::mutate(id_comuna = sprintf("%05d", id_comuna))

  # Final type casting
  all_series <- all_series |>
    dplyr::mutate(id_comuna = as.character(id_comuna)) |>
    dplyr::mutate(
      col01 = suppressWarnings(as.integer(col01)),
      col02 = suppressWarnings(as.integer(col02)),
      col03 = suppressWarnings(as.integer(col03)),
      col04 = suppressWarnings(as.integer(col04))
    ) |>
    dplyr::mutate_if(is.numeric, as.integer) |>
    dplyr::mutate_if(is.logical, as.integer)

  return(all_series)
}


# -------------------------------------------------------------------------
# make_all_series_diabetes()
# Joins DPA lookup, corrects mismatched id_comunas, filters for diabetes
# codes and December only.
# -------------------------------------------------------------------------

make_all_series_diabetes <- function(all_series_data, dpa_data) {

  # Step 1: Build mismatch lookup (old id_comuna → new id_comuna)
  dpa_mismatches_renamed <- dpa_data |>
    dplyr::filter(id_comuna != id_comuna2) |>
    dplyr::select(id_comuna_old = id_comuna, id_comuna_new = id_comuna2)

  # Step 2: Replace mismatched id_comunas
  all_series_updated <- all_series_data |>
    dplyr::left_join(dpa_mismatches_renamed, by = c("id_comuna" = "id_comuna_old")) |>
    dplyr::mutate(id_comuna = dplyr::coalesce(id_comuna_new, id_comuna)) |>
    dplyr::select(-id_comuna_new)

  # Step 3: Build final comuna lookup from DPA (keyed on updated id_comuna2)
  comuna_lookup <- dpa_data |>
    dplyr::select(id_comuna = id_comuna2, comuna, region2, id_region2, servicio_salud, id_servicio2) |>
    dplyr::distinct()

  # Step 4: Join, update region/service columns, filter for diabetes codes + December
  all_series_diabetes <- all_series_updated |>
    dplyr::left_join(comuna_lookup, by = "id_comuna") |>
    dplyr::mutate(
      id_region  = id_region2,
      region     = region2,
      id_servicio = id_servicio2
    ) |>
    dplyr::filter(
      codigo_prestacion %in% c("P4150602", "P4190950", "P4190400", "P4180300"),
      mes == 12
    ) |>
    dplyr::select(-region2, -id_region2, -id_servicio2)

  return(all_series_diabetes)
}


# -------------------------------------------------------------------------
# make_coverage_data()
# Aggregates to comuna level, computes coverage ratios, removes Q1 comunas.
# -------------------------------------------------------------------------

make_coverage_data <- function(all_series_diabetes) {

  # Add geographic zone variable
  coverage <- all_series_diabetes |>
    dplyr::mutate(zona = dplyr::case_when(
      region %in% c("De Arica y Parinacota", "De Tarapacá",
                    "De Antofagasta", "De Atacama", "De Coquimbo") ~ 1,
      region %in% c("De Valparaíso", "Metropolitana de Santiago",
                    "Del Libertador B. O'Higgins", "Del Maule")     ~ 2,
      TRUE                                                           ~ 3
    ))

  # Aggregate to comuna × year × codigo_prestacion, then pivot wide
  coverage <- coverage |>
    dplyr::group_by(ano, comuna, id_comuna, id_region, region, zona, codigo_prestacion) |>
    dplyr::summarise(cantidad = round(sum(col01)), .groups = "drop") |>
    tidyr::spread(codigo_prestacion, cantidad) |>
    dplyr::rename(
      dm            = P4150602,
      dm_fo         = P4190950,
      dm_fo_2       = P4190400,
      dm_hg_menor7  = P4180300
    ) |>
    # Coalesce the two diabetic retinopathy screening codes
    dplyr::mutate(dm_fo = dplyr::coalesce(dm_fo, dm_fo_2)) |>
    dplyr::select(-dm_fo_2) |>
    dplyr::mutate(
      drs_coverage = dm_fo / dm,
      dm_coverage  = dm_hg_menor7 / dm
    )

  # Remove Q1 comunas (smallest diabetes burden — pooled across all years)
  # NOTE: quintiles are computed on the full pooled dataset, so a municipality
  # is excluded entirely if it falls in Q1 in the pooled distribution.
  coverage_noq1 <- coverage |>
    dplyr::ungroup() |>
    dplyr::mutate(quintil_dm_category = cut(
      dm,
      breaks = quantile(dm, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE),
      labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
      include.lowest = TRUE
    )) |>
    dplyr::filter(quintil_dm_category != "Q1") |>
    dplyr::mutate(
      comuna2      = match(comuna, unique(comuna)),
      ano2         = ano - 2011,
      drs_coverage = replace(drs_coverage, drs_coverage > 1, 1)
    ) |>
    dplyr::arrange(ano, comuna)

  # Final modelling dataset
  coverage_data <- coverage_noq1 |>
    dplyr::arrange(comuna2, ano) |>
    dplyr::mutate(
      year = ano - min(ano),
      id   = comuna2,
      drsc = drs_coverage,
      dgcc = dm_coverage
    )

  return(coverage_data)
}


# -------------------------------------------------------------------------
# make_model_adequacy_table()
# Builds the full model comparison / adequacy table from all_named_models.
# -------------------------------------------------------------------------

make_model_adequacy_table <- function(all_named_models) {

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
      tk <- LCTMtoolkit::LCTMtoolkit(model)
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
        Model                       = model_name,
        Smallest_Class_Size_Percentage = pp$size,
        Smallest_Class_Count        = pp$count,
        Lowest_OCC                  = occ$occ,
        Lowest_APPA                 = occ$appa,
        Highest_Mismatch            = occ$mismatch,
        VLMRLRT_P_Value             = p_val
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

  # --- Helper: Mahalanobis degree of separation ---------------------------

  extract_coeffs_by_class <- function(model) {
    tryCatch({
      if (model$ng <= 1) return(NULL)
      beta <- model$best
      df   <- data.frame(
        intercept = beta[1:model$ng],
        slope     = beta[(model$ng + 1):(2 * model$ng)]
      )
      if (any(is.na(df))) return(NULL)
      df
    }, error = function(e) NULL)
  }

  calculate_Mahalanobis_DoS <- function(df, cov_mat = NULL) {
    if (is.null(df) || nrow(df) <= 1) return(NA_real_)
    mat     <- as.matrix(df)
    cov_mat <- if (is.null(cov_mat)) diag(ncol(mat)) else cov_mat
    inv_cov <- tryCatch(solve(cov_mat), error = function(e) NULL)
    if (is.null(inv_cov)) return(NA_real_)
    dist <- combn(1:nrow(mat), 2, function(ix) {
      d <- mat[ix[1], ] - mat[ix[2], ]
      sqrt(t(d) %*% inv_cov %*% d)
    })
    mean(dist)
  }

  add_Mahalanobis_DoS <- function(model_list, table) {
    df <- lapply(names(model_list), function(name) {
      coeff <- extract_coeffs_by_class(model_list[[name]])
      dos   <- calculate_Mahalanobis_DoS(coeff)
      tibble::tibble(Model = name, DoS_Mahalanobis = round(dos, 4))
    }) |> dplyr::bind_rows()
    dplyr::left_join(table, df, by = "Model")
  }

  # --- Run the pipeline ---------------------------------------------------

  summary_table   <- build_summary_table(names(all_named_models), all_named_models)
  diagnostic_table <- process_all_models(all_named_models)

  summary_table$Model    <- as.character(summary_table$Model)
  diagnostic_table$Model <- as.character(diagnostic_table$Model)

  adequacy_table <- dplyr::left_join(summary_table, diagnostic_table, by = "Model")
  final_table    <- add_Mahalanobis_DoS(all_named_models, adequacy_table)

  # Tag model structure
  class_cols <- grep("^%class\\d+$", names(final_table), value = TRUE)

  model_adequacy_table <- final_table |>
    dplyr::mutate(structure = dplyr::case_when(
      stringr::str_detect(Model, "linear_nre_homocedastic")    ~ "A",
      stringr::str_detect(Model, "linear_nre_heterocedastic")  ~ "B",
      stringr::str_detect(Model, "quadratic_nre")              ~ "C",
      stringr::str_detect(Model, "cubic_nre")                  ~ "D",
      stringr::str_detect(Model, "linear_random_intercept_slope") ~ "F",
      stringr::str_detect(Model, "linear_random_intercept")    ~ "E",
      stringr::str_detect(Model, "quadratic_random_effects_prop") ~ "H",
      stringr::str_detect(Model, "quadratic_random_effects")   ~ "G",
      stringr::str_detect(Model, "cubic_random_effects_prop")  ~ "J",
      stringr::str_detect(Model, "cubic_random_effects")       ~ "I",
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
