# R/00_data.R
# -------------------------------------------------------------------------
# Data construction pipeline for municipal diabetes coverage analysis
#   - make_all_series():           combine annual SerieP files (2011-2023)
#   - make_all_series_diabetes():  join DPA lookup, filter diabetes codes
#                                  and December observations
#   - make_coverage_data():        aggregate to comuna level, compute
#                                  coverage ratios, apply exclusion criteria
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# make_all_series()
# Combines a list of yearly data frames and applies final cleaning/type
# casting. All inputs must be character-only for safe row-binding across
# years (years have heterogeneous column types upstream).
# -------------------------------------------------------------------------

make_all_series <- function(series_list) {
  
  all_series <- purrr::list_rbind(series_list)
  
  all_series <- type.convert(all_series, as.is = TRUE)
  
  all_series <- all_series |>
    dplyr::mutate(id_comuna = sprintf("%05d", id_comuna))
  
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
# Joins DPA lookup table, corrects pre-2018 mismatched id_comunas, filters
# for the four diabetes-related prestacion codes and December observations
# only. December captures cumulative annual registrations under PSCV.
# -------------------------------------------------------------------------

make_all_series_diabetes <- function(all_series_data, dpa_data) {
  
  # Build mismatch lookup (old id_comuna → new id_comuna)
  dpa_mismatches_renamed <- dpa_data |>
    dplyr::filter(id_comuna != id_comuna2) |>
    dplyr::select(id_comuna_old = id_comuna, id_comuna_new = id_comuna2)
  
  # Replace mismatched id_comunas
  all_series_updated <- all_series_data |>
    dplyr::left_join(dpa_mismatches_renamed, by = c("id_comuna" = "id_comuna_old")) |>
    dplyr::mutate(id_comuna = dplyr::coalesce(id_comuna_new, id_comuna)) |>
    dplyr::select(-id_comuna_new)
  
  # Build final comuna lookup keyed on updated id_comuna2
  comuna_lookup <- dpa_data |>
    dplyr::select(id_comuna = id_comuna2, comuna, region2, id_region2,
                  servicio_salud, id_servicio2) |>
    dplyr::distinct()
  
  # Join, update region/service columns, filter for diabetes codes + December
  all_series_diabetes <- all_series_updated |>
    dplyr::left_join(comuna_lookup, by = "id_comuna") |>
    dplyr::mutate(
      id_region   = id_region2,
      region      = region2,
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
# Aggregates raw municipal-level diabetes records to municipality × year,
# computes coverage ratios (EGCC and DRSC), and applies two exclusion
# criteria for trajectory estimation:
#
#   1. Pooled-Q1 filter: removes municipality-year observations in the
#      lowest quintile of annual diabetes registrations to reduce
#      instability from sparse denominators.
#
#   2. Structural filter: removes any municipality with fewer than two
#      remaining annual observations after the Q1 filter, because
#      trajectory estimation requires at least two repeated measurements.
#
# Yields 299 analytical municipalities (the manuscript denominator).
#
# Note on variable naming: internally `dgcc` corresponds to the manuscript
# indicator `EGCC` (effective glycaemic control coverage); `dm_coverage`
# is its raw computation. `drsc` corresponds to `DRSC`.
# -------------------------------------------------------------------------

make_coverage_data <- function(all_series_diabetes) {
  
  coverage <- all_series_diabetes |>
    dplyr::mutate(zona = dplyr::case_when(
      region %in% c("De Arica y Parinacota", "De Tarapacá",
                    "De Antofagasta", "De Atacama", "De Coquimbo") ~ 1,
      region %in% c("De Valparaíso", "Metropolitana de Santiago",
                    "Del Libertador B. O'Higgins", "Del Maule")     ~ 2,
      TRUE                                                           ~ 3
    ))
  
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
    dplyr::mutate(dm_fo = dplyr::coalesce(dm_fo, dm_fo_2)) |>
    dplyr::select(-dm_fo_2) |>
    dplyr::mutate(
      drs_coverage = dm_fo / dm,
      dm_coverage  = dm_hg_menor7 / dm
    )
  
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
  
  coverage_noq1 <- coverage_noq1 |>
    dplyr::group_by(comuna2) |>
    dplyr::filter(dplyr::n() >= 2) |>
    dplyr::ungroup()
  
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