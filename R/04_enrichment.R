# R/04_enrichment.R
# -------------------------------------------------------------------------
# Enrichment of coverage_data with municipal-level covariates and latent
# class assignments. Produces the analytic data frame used downstream in
# multinomial regression, mixed-effects models, and cross-classification.
#
# Functions:
#   - norm_comuna():                       commune name normaliser
#   - make_coverage_data_enriched():       coverage_data + IDC + PNDR
#                                           + zone (Northern/Centre/Southern)
#                                           + urbanicity (Rural/Mixed/Urban)
#                                           + Melipeuco exclusion
#   - make_coverage_classes():             extract posterior class assignments
#                                           and probabilities for both
#                                           outcomes from fitted models
#   - make_coverage_summary():             municipality-level summary stats
#                                           (means, SDs, change metrics)
#   - make_coverage_columns_wide():        wide-format coverage by year
#                                           (used downstream for tables)
#   - make_coverage_demographics_class():  final analytic table joining
#                                           summary + demographics + classes
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# norm_comuna()
# Normalises commune names: trim, squish, uppercase, strip diacritics.
# -------------------------------------------------------------------------

norm_comuna <- function(x) {
  x |>
    stringr::str_trim() |>
    stringr::str_squish() |>
    stringr::str_to_upper() |>
    stringi::stri_trans_general("Latin-ASCII")
}


# -------------------------------------------------------------------------
# make_coverage_data_enriched()
# Joins IDC and PNDR to coverage_data using normalised commune keys and
# a small alias map (AISEN/COIHAIQUE/CALERA/NATALES). Recodes zone and
# urbanicity to manuscript labels. Excludes Melipeuco.
# -------------------------------------------------------------------------

make_coverage_data_enriched <- function(coverage_data, idc, comunas_pndr) {
  
  coverage_clean <- coverage_data |>
    dplyr::mutate(comuna_key = norm_comuna(comuna))
  
  idc_clean <- idc |>
    dplyr::mutate(comuna_key = norm_comuna(comuna))
  
  idc_keys <- idc_clean |>
    dplyr::distinct(comuna_key) |>
    dplyr::pull(comuna_key)
  
  alias_map <- tibble::tibble(
    from = c("AISEN", "COIHAIQUE", "CALERA", "NATALES"),
    to   = c(
      if ("AYSEN"          %in% idc_keys) "AYSEN"          else "AISEN",
      if ("COYHAIQUE"      %in% idc_keys) "COYHAIQUE"      else "COIHAIQUE",
      if ("LA CALERA"      %in% idc_keys) "LA CALERA"      else "CALERA",
      if ("PUERTO NATALES" %in% idc_keys) "PUERTO NATALES" else "NATALES"
    )
  )
  
  coverage_fix <- coverage_clean |>
    dplyr::left_join(alias_map, by = c("comuna_key" = "from")) |>
    dplyr::mutate(comuna_key_fix = dplyr::coalesce(to, comuna_key)) |>
    dplyr::select(-to)
  
  coverage_idc <- coverage_fix |>
    dplyr::left_join(
      idc_clean |> dplyr::select(comuna_key, region:rangos, comuna_idc = comuna),
      by = c("comuna_key_fix" = "comuna_key"),
      suffix = c("_cov", "_idc")
    )
  
  main <- coverage_idc |>
    dplyr::mutate(join_key = dplyr::coalesce(comuna_key_fix, comuna_key))
  
  target_keys <- unique(main$join_key)
  
  pndr_clean <- comunas_pndr |>
    dplyr::mutate(join_key = norm_comuna(comuna)) |>
    dplyr::mutate(
      join_key = dplyr::case_when(
        join_key == "CALERA"    & "LA CALERA"      %in% target_keys ~ "LA CALERA",
        join_key == "NATALES"   & "PUERTO NATALES" %in% target_keys ~ "PUERTO NATALES",
        join_key == "COIHAIQUE" & "COYHAIQUE"      %in% target_keys ~ "COYHAIQUE",
        TRUE ~ join_key
      )
    ) |>
    dplyr::distinct(join_key, .keep_all = TRUE) |>
    dplyr::rename(comuna_pndr = comuna)
  
  coverage_idc_pndr <- main |>
    dplyr::left_join(
      pndr_clean |>
        dplyr::select(
          join_key, comuna_pndr,
          cod_reg, region_pndr = region, cod_com,
          n_habitantes, km_2, densidad, tipo_com, clasificacion
        ),
      by = "join_key"
    )
  
  coverage_data_enriched <- coverage_idc_pndr |>
    dplyr::mutate(
      zone = dplyr::case_when(
        zona == 1 ~ "Northern",
        zona == 2 ~ "Centre",
        zona == 3 ~ "Southern",
        TRUE      ~ NA_character_
      ),
      urbanicity = dplyr::case_when(
        clasificacion == "Urbana" ~ "Urban",
        clasificacion == "Mixta"  ~ "Mixed",
        clasificacion == "Rural"  ~ "Rural",
        TRUE                      ~ NA_character_
      )
    ) |>
    dplyr::filter(comuna != "Melipeuco")
  
  return(coverage_data_enriched)
}


# -------------------------------------------------------------------------
# make_coverage_classes()
# Extracts posterior class assignments and probabilities from the two
# retained 4-class cubic NRE models. Adds substantive labels for EGCC
# and DRSC trajectory classes.
# -------------------------------------------------------------------------

make_coverage_classes <- function(all_named_models) {
  
  dgcc_class <- all_named_models[["4class_cubic_nre_dgcc_model"]]$pprob |>
    dplyr::rename(
      class_dgcc = class,
      prob1_dgcc = prob1,
      prob2_dgcc = prob2,
      prob3_dgcc = prob3,
      prob4_dgcc = prob4
    )
  
  drsc_class <- all_named_models[["4class_cubic_nre_drsc_model"]]$pprob |>
    dplyr::rename(
      class_drsc = class,
      prob1_drsc = prob1,
      prob2_drsc = prob2,
      prob3_drsc = prob3,
      prob4_drsc = prob4
    )
  
  coverage_classes <- dgcc_class |>
    dplyr::left_join(drsc_class, by = "id") |>
    dplyr::mutate(
      class_dgcc_label = dplyr::case_when(
        class_dgcc == 1 ~ "Lowest",
        class_dgcc == 2 ~ "Highest",
        class_dgcc == 3 ~ "Stable upper-medium",
        class_dgcc == 4 ~ "Stable medium",
        TRUE            ~ NA_character_
      ),
      class_drsc_label = dplyr::case_when(
        class_drsc == 1 ~ "Increasing",
        class_drsc == 2 ~ "Stable medium",
        class_drsc == 3 ~ "Predominantly highest",
        class_drsc == 4 ~ "Lowest",
        TRUE            ~ NA_character_
      )
    )
  
  return(coverage_classes)
}


# -------------------------------------------------------------------------
# make_coverage_summary()
# Per-municipality summary of coverage variables: counts of valid years,
# means, SDs, medians, first/last non-missing values, absolute and
# relative change over the period.
# -------------------------------------------------------------------------

make_coverage_summary <- function(coverage_data_enriched) {
  
  coverage_data_enriched |>
    dplyr::group_by(comuna) |>
    dplyr::arrange(ano, .by_group = TRUE) |>
    dplyr::summarise(
      n_years   = dplyr::n_distinct(ano),
      n_dm      = sum(!is.na(dm)),
      n_drsc    = sum(!is.na(drsc)),
      n_dgcc    = sum(!is.na(dgcc)),
      n_all3    = sum(!is.na(dm) & !is.na(drsc) & !is.na(dgcc)),
      prop_dm   = n_dm   / n_years,
      prop_drsc = n_drsc / n_years,
      prop_dgcc = n_dgcc / n_years,
      mean_dm   = mean(dm,   na.rm = TRUE),
      sd_dm     = sd(dm,     na.rm = TRUE),
      median_drsc = median(drsc, na.rm = TRUE),
      mean_drsc   = mean(drsc,   na.rm = TRUE),
      sd_drsc     = sd(drsc,     na.rm = TRUE),
      median_dgcc = median(dgcc, na.rm = TRUE),
      mean_dgcc   = mean(dgcc,   na.rm = TRUE),
      sd_dgcc     = sd(dgcc,     na.rm = TRUE),
      first_drsc      = dplyr::first(stats::na.omit(drsc)),
      last_drsc       = dplyr::last(stats::na.omit(drsc)),
      change_drsc     = last_drsc - first_drsc,
      rel_change_drsc = dplyr::if_else(
        is.na(first_drsc) | first_drsc == 0, NA_real_,
        100 * change_drsc / first_drsc
      ),
      first_dgcc      = dplyr::first(stats::na.omit(dgcc)),
      last_dgcc       = dplyr::last(stats::na.omit(dgcc)),
      change_dgcc     = last_dgcc - first_dgcc,
      rel_change_dgcc = dplyr::if_else(
        is.na(first_dgcc) | first_dgcc == 0, NA_real_,
        100 * change_dgcc / first_dgcc
      ),
      .groups = "drop"
    )
}


# -------------------------------------------------------------------------
# make_coverage_columns_wide()
# Wide-format coverage by year for a given outcome variable. Used as
# input for descriptive tables outside the pipeline.
# -------------------------------------------------------------------------

make_coverage_columns_wide <- function(coverage_data_enriched, coverage_var) {
  
  cov_name <- rlang::as_name(rlang::ensym(coverage_var))
  
  coverage_data_enriched |>
    dplyr::group_by(comuna, year) |>
    dplyr::summarise(
      mean = mean({{ coverage_var }}, na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from  = year,
      values_from = mean,
      names_glue  = "{.value}_{year}"
    ) |>
    dplyr::rename_with(~ paste0(cov_name, "_", .x), -comuna)
}


# -------------------------------------------------------------------------
# make_coverage_demographics_class()
# Final analytic data frame combining: per-municipality summary stats,
# demographic covariates (IDC + PNDR), and latent class assignments
# with substantive labels. This is the input for downstream multinomial
# regression and cross-classification analyses.
# -------------------------------------------------------------------------

make_coverage_demographics_class <- function(coverage_data_enriched,
                                             coverage_classes) {
  
  coverage_summary <- make_coverage_summary(coverage_data_enriched)
  
  additional_data <- coverage_data_enriched |>
    dplyr::select(
      comuna, id, id_comuna, id_region, region_cov,
      zone, urbanicity,
      bienestar, economia, educacion, idc, ranking,
      n_habitantes, km_2, densidad
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(scaled_idc = as.numeric(scale(idc)))
  
  coverage_summary_demographics <- coverage_summary |>
    dplyr::left_join(additional_data, by = "comuna")
  
  coverage_summary_demographics |>
    dplyr::left_join(coverage_classes, by = "id")
}