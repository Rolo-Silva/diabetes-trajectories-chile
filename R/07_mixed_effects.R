# R/07_mixed_effects.R
# -------------------------------------------------------------------------
# Ecological linear mixed-effects analysis of EGCC and DRSC trajectories
# across Chilean municipalities (2011–2023).
#
# Two model specifications per outcome:
#   - MAIN:        outcome ~ year_c + idc_mean_z + urbanicity + zone + (1 | id)
#                   year_c = average annual change shared across municipalities
#                   covariates = differences in mean level
#
#   - INTERACTION: outcome ~ year_c * (idc_mean_z + urbanicity + zone) + (1 | id)
#                   year_c x covariate = differences in annual change (slope)
#
# Models use random intercept by municipality. Time is centred at 2011
# (year_c = ano - 2011). IDC is standardised at the municipality level
# (mean IDC per municipality, then z-scored across municipalities) so
# the coefficient represents the difference per 1 SD of deprivation.
#
# Output: list with the four model fits, four tidy coefficient tables
# with Wald 95% CIs, two LRT comparisons (main vs interaction), and a
# preformatted Supplementary Table 7 (interaction-model-only, both
# outcomes stacked).
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# .me_extract_lmer_results()
# Internal: extract fixed-effects coefficients + Wald 95% CIs + p-values
# from an lmer fit. Returns a tidy data frame; avoids broom.mixed::tidy
# inconsistencies across versions.
# -------------------------------------------------------------------------

.me_extract_lmer_results <- function(model) {
  
  coefs       <- as.data.frame(summary(model)$coefficients)
  coefs$term  <- rownames(coefs)
  rownames(coefs) <- NULL
  
  # Defensive: handle both lmerMod (no p-values) and lmerModLmerTest
  if (!"Pr(>|t|)" %in% names(coefs)) {
    coefs[["Pr(>|t|)"]] <- NA_real_
  }
  if (!"t value" %in% names(coefs)) {
    coefs[["t value"]] <- coefs$Estimate / coefs[["Std. Error"]]
  }
  
  ci <- as.data.frame(stats::confint(model, parm = "beta_", method = "Wald"))
  ci$term     <- rownames(ci)
  rownames(ci) <- NULL
  names(ci)[1:2] <- c("conf.low", "conf.high")
  
  coefs |>
    dplyr::rename(
      estimate  = Estimate,
      std.error = `Std. Error`,
      statistic = `t value`,
      p.value   = `Pr(>|t|)`
    ) |>
    dplyr::left_join(ci, by = "term") |>
    dplyr::select(term, estimate, std.error, statistic, p.value,
                  conf.low, conf.high)
}


# -------------------------------------------------------------------------
# .me_build_long_data()
# Internal: convert coverage_data_enriched into a municipality-year long
# panel ready for lmer, with centred year, harmonised factor levels for
# zone and urbanicity, and municipality-level standardised IDC.
# -------------------------------------------------------------------------

.me_build_long_data <- function(coverage_data_enriched) {
  
  df <- coverage_data_enriched |>
    dplyr::transmute(
      id         = as.numeric(id),
      id_comuna  = as.character(id_comuna),
      comuna     = as.character(comuna),
      ano        = as.numeric(ano),
      dgcc       = as.numeric(dgcc),
      drsc       = as.numeric(drsc),
      idc        = as.numeric(idc),
      urbanicity = dplyr::case_when(
        urbanicity == "Rural" ~ "Rural",
        urbanicity == "Mixed" ~ "Mixed",
        urbanicity == "Urban" ~ "Urban",
        TRUE ~ NA_character_
      ),
      zone = dplyr::case_when(
        zone == "Northern" ~ "Northern",
        zone == "Centre"   ~ "Centre",
        zone == "Southern" ~ "Southern",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::mutate(
      urbanicity = factor(urbanicity, levels = c("Mixed", "Rural", "Urban")),
      zone       = factor(zone,       levels = c("Centre", "Northern", "Southern")),
      year_c     = ano - min(ano, na.rm = TRUE)
    ) |>
    dplyr::arrange(id, ano)
  
  # Municipality-level standardised IDC
  idc_by_commune <- df |>
    dplyr::group_by(id) |>
    dplyr::summarise(idc_mean = mean(idc, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(idc_mean_z = as.numeric(scale(idc_mean)))
  
  df |>
    dplyr::select(-idc) |>
    dplyr::left_join(idc_by_commune, by = "id")
}


# -------------------------------------------------------------------------
# .me_format_supp_table_row()
# Internal: convert tidy lmer results into ST7 row format.
# -------------------------------------------------------------------------

.me_format_supp_table_row <- function(res_df, outcome_label, terms_keep,
                                      term_labels) {
  
  res_df |>
    dplyr::filter(term %in% terms_keep) |>
    dplyr::mutate(
      Characteristic = term_labels[term],
      beta_ci = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high),
      p_value = dplyr::if_else(
        p.value < 0.001, "<0.001",
        sprintf("%.3f", p.value)
      ),
      Outcome = outcome_label
    ) |>
    dplyr::select(Outcome, Characteristic, beta_ci, p_value)
}


# -------------------------------------------------------------------------
# make_mixed_effects_analysis()
# Main orchestrator. Fits MAIN and INTERACTION models for both outcomes,
# extracts tidy coefficient tables, runs LRT comparisons, and produces
# preformatted Supplementary Table 7.
# -------------------------------------------------------------------------

make_mixed_effects_analysis <- function(coverage_data_enriched) {
  
  df_long <- .me_build_long_data(coverage_data_enriched)
  
  # ---- Fit four models ----
  
  mod_dgcc_main <- lmerTest::lmer(
    dgcc ~ year_c + idc_mean_z + urbanicity + zone + (1 | id),
    data = df_long, REML = FALSE
  )
  
  mod_drsc_main <- lmerTest::lmer(
    drsc ~ year_c + idc_mean_z + urbanicity + zone + (1 | id),
    data = df_long, REML = FALSE
  )
  
  mod_dgcc_int <- lmerTest::lmer(
    dgcc ~ year_c * (idc_mean_z + urbanicity + zone) + (1 | id),
    data = df_long, REML = FALSE
  )
  
  mod_drsc_int <- lmerTest::lmer(
    drsc ~ year_c * (idc_mean_z + urbanicity + zone) + (1 | id),
    data = df_long, REML = FALSE
  )
  
  # ---- Extract tidy results ----
  
  res_dgcc_main <- .me_extract_lmer_results(mod_dgcc_main)
  res_drsc_main <- .me_extract_lmer_results(mod_drsc_main)
  res_dgcc_int  <- .me_extract_lmer_results(mod_dgcc_int)
  res_drsc_int  <- .me_extract_lmer_results(mod_drsc_int)
  
  # ---- LRT comparisons (main vs interaction) ----
  
  lrt_dgcc <- stats::anova(mod_dgcc_main, mod_dgcc_int)
  lrt_drsc <- stats::anova(mod_drsc_main, mod_drsc_int)
  
  lrt_notes <- tibble::tibble(
    Outcome = c("EGCC", "DRSC"),
    LRT_chisq = c(lrt_dgcc$Chisq[2], lrt_drsc$Chisq[2]),
    LRT_df    = c(lrt_dgcc$Df[2],    lrt_drsc$Df[2]),
    LRT_p     = c(lrt_dgcc$`Pr(>Chisq)`[2], lrt_drsc$`Pr(>Chisq)`[2]),
    LRT_p_fmt = dplyr::if_else(
      LRT_p < 0.001, "<0.001",
      sprintf("%.3f", LRT_p)
    )
  )
  
  # ---- Build Supplementary Table 7 (interaction-only) ----
  
  interaction_terms <- c(
    "year_c",
    "idc_mean_z",
    "urbanicityRural",
    "urbanicityUrban",
    "zoneNorthern",
    "zoneSouthern",
    "year_c:idc_mean_z",
    "year_c:urbanicityRural",
    "year_c:urbanicityUrban",
    "year_c:zoneNorthern",
    "year_c:zoneSouthern"
  )
  
  interaction_labels <- c(
    "year_c"                  = "Year (annual change)",
    "idc_mean_z"              = "ICD (per 1 SD higher)",
    "urbanicityRural"         = "Rural vs mixed",
    "urbanicityUrban"         = "Urban vs mixed",
    "zoneNorthern"            = "Northern vs centre",
    "zoneSouthern"            = "Southern vs centre",
    "year_c:idc_mean_z"       = "Year x ICD",
    "year_c:urbanicityRural"  = "Year x rural",
    "year_c:urbanicityUrban"  = "Year x urban",
    "year_c:zoneNorthern"     = "Year x northern",
    "year_c:zoneSouthern"     = "Year x southern"
  )
  
  supp_table7 <- dplyr::bind_rows(
    .me_format_supp_table_row(res_dgcc_int, "EGCC",
                              interaction_terms, interaction_labels),
    .me_format_supp_table_row(res_drsc_int, "DRSC",
                              interaction_terms, interaction_labels)
  )
  
  # ---- Return all components ----
  
  list(
    df_long       = df_long,
    models = list(
      dgcc_main = mod_dgcc_main,
      drsc_main = mod_drsc_main,
      dgcc_int  = mod_dgcc_int,
      drsc_int  = mod_drsc_int
    ),
    results = list(
      dgcc_main = res_dgcc_main,
      drsc_main = res_drsc_main,
      dgcc_int  = res_dgcc_int,
      drsc_int  = res_drsc_int
    ),
    lrt        = list(dgcc = lrt_dgcc, drsc = lrt_drsc),
    lrt_notes  = lrt_notes,
    supp_table7 = supp_table7
  )
}

