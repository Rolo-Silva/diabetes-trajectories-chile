# R/05_multinom.R
# -------------------------------------------------------------------------
# Multinomial logistic regression with pseudo-class draws and Rubin's
# rules pooling for trajectory class memberships.
#
# For each outcome (DGCC, DRSC) and each reference class (1..4), fits
# a multinomial model with predictors zone + urbanicity + scaled_idc.
# Uses posterior class probabilities to perform multiple imputation of
# class assignment, then pools coefficients via Rubin's rules.
#
# Functions:
#   - run_multinom_MI():           single multinomial fit with MI pooling
#                                   for one outcome and one reference class
#   - make_multinom_validation():  orchestrator: runs all 4 ref_class × 2
#                                   outcomes = 8 tables and returns them
#                                   as a named list ready for downstream
#                                   reporting
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# run_multinom_MI()
# Fits multinomial logistic regression with pseudo-class multiple
# imputation. Configurable reference classes for outcome, zone, and
# urbanicity. Returns a wide-format tibble of RRR (95% CI).
#
# Arguments:
#   data:         data frame with class probabilities + covariates
#   prob_cols:    character vector of 4 posterior probability columns
#                  e.g. c("prob1_dgcc","prob2_dgcc","prob3_dgcc","prob4_dgcc")
#   outcome_name: "DGCC" or "DRSC" (for attribute tagging only)
#   n_imp:        number of imputed datasets (default 10)
#   ref_class:    reference class for outcome (1..4)
#   zone_ref:     reference level for zone (e.g. "Centre"). NULL = first level.
#   urban_ref:    reference level for urbanicity (e.g. "Mixed"). NULL = first.
# -------------------------------------------------------------------------

run_multinom_MI <- function(data,
                            prob_cols,
                            outcome_name,
                            n_imp     = 10,
                            ref_class = 1,
                            zone_ref  = NULL,
                            urban_ref = NULL) {
  
  n_classes <- length(prob_cols)
  stopifnot(ref_class %in% seq_len(n_classes))
  stopifnot(all(prob_cols %in% names(data)))
  
  set.seed(123)
  
  # ---- 1) Pseudo-class multiple imputation ----
  
  imputed_datasets <- lapply(seq_len(n_imp), function(i) {
    
    prob_mat <- as.matrix(data[, prob_cols])
    data$class_imp <- apply(prob_mat, 1, function(p)
      sample(seq_len(n_classes), 1, prob = p)
    )
    
    # Reference for outcome
    levs <- as.character(c(ref_class, setdiff(seq_len(n_classes), ref_class)))
    data$class_fac <- factor(data$class_imp, levels = levs)
    
    # Reference for predictors (if specified and present)
    if (!is.null(zone_ref) && "zone" %in% names(data)) {
      data$zone <- factor(
        data$zone,
        levels = c(zone_ref, setdiff(unique(stats::na.omit(data$zone)), zone_ref))
      )
    }
    if (!is.null(urban_ref) && "urbanicity" %in% names(data)) {
      data$urbanicity <- factor(
        data$urbanicity,
        levels = c(urban_ref, setdiff(unique(stats::na.omit(data$urbanicity)), urban_ref))
      )
    }
    
    data
  })
  
  # ---- 2) Fit multinomial logistic regression on each imputed dataset ----
  
  models <- lapply(imputed_datasets, function(d) {
    nnet::multinom(
      class_fac ~ zone + urbanicity + scaled_idc,
      data  = d,
      trace = FALSE
    )
  })
  
  # ---- 3) Tidy each model ----
  
  tidy_list <- lapply(models, broom::tidy)
  
  # ---- 4) Pool coefficients with Rubin's rules ----
  
  pool_rubins_rules <- function(tidy_list) {
    all <- dplyr::bind_rows(tidy_list, .id = "imp")
    all |>
      dplyr::group_by(y.level, term) |>
      dplyr::summarise(
        mean_est = mean(estimate, na.rm = TRUE),
        mean_se  = sqrt(
          mean(std.error^2, na.rm = TRUE) +
            (1 + 1 / length(tidy_list)) * stats::var(estimate, na.rm = TRUE)
        ),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        RRR   = exp(mean_est),
        lower = exp(mean_est - 1.96 * mean_se),
        upper = exp(mean_est + 1.96 * mean_se)
      )
  }
  
  pooled <- pool_rubins_rules(tidy_list)
  
  # ---- 5) Relabel terms and pivot to wide ----
  
  table_out <- pooled |>
    dplyr::mutate(
      RRR_CI = sprintf("%.2f (%.2f–%.2f)", RRR, lower, upper),
      term = dplyr::case_when(
        term == "(Intercept)" ~ "Intercept",
        term == "scaled_idc"  ~ "Index of Deprivation (scaled)",
        term == "urbanicityRural" & !is.null(urban_ref) & urban_ref == "Mixed" ~ "Urbanicity: Rural vs Mixed",
        term == "urbanicityUrban" & !is.null(urban_ref) & urban_ref == "Mixed" ~ "Urbanicity: Urban vs Mixed",
        term == "zoneNorthern"   & !is.null(zone_ref)  & zone_ref == "Centre"  ~ "Zone: Northern vs Centre",
        term == "zoneSouthern"   & !is.null(zone_ref)  & zone_ref == "Centre"  ~ "Zone: Southern vs Centre",
        TRUE ~ term
      ),
      comp_col = paste0("Class ", y.level, " vs Class ", ref_class)
    ) |>
    dplyr::select(term, comp_col, RRR_CI) |>
    tidyr::pivot_wider(
      names_from  = comp_col,
      values_from = RRR_CI,
      values_fill = NA
    )
  
  attr(table_out, "outcome")   <- outcome_name
  attr(table_out, "ref_class") <- ref_class
  
  table_out
}


# -------------------------------------------------------------------------
# make_multinom_validation()
# Orchestrator: runs all 4 reference classes × 2 outcomes = 8 multinomial
# tables. Returns a named list with all tables.
#
# Naming convention: "<outcome>_ref<class>", e.g. "dgcc_ref1", "drsc_ref4".
#
# Arguments:
#   coverage_demographics_class: analytic data frame with class probs
#                                  and covariates (zone, urbanicity, scaled_idc)
#   zone_ref:                    reference level for zone (default "Centre")
#   urban_ref:                   reference level for urbanicity (default "Mixed")
#   n_imp:                       imputations per fit (default 10)
# -------------------------------------------------------------------------

make_multinom_validation <- function(coverage_demographics_class,
                                     zone_ref  = "Centre",
                                     urban_ref = "Mixed",
                                     n_imp     = 10) {
  
  prob_cols_dgcc <- c("prob1_dgcc", "prob2_dgcc", "prob3_dgcc", "prob4_dgcc")
  prob_cols_drsc <- c("prob1_drsc", "prob2_drsc", "prob3_drsc", "prob4_drsc")
  
  ref_classes <- 1:4
  
  tables_dgcc <- purrr::map(ref_classes, function(rc) {
    run_multinom_MI(
      data         = coverage_demographics_class,
      prob_cols    = prob_cols_dgcc,
      outcome_name = "DGCC",
      n_imp        = n_imp,
      ref_class    = rc,
      zone_ref     = zone_ref,
      urban_ref    = urban_ref
    )
  })
  names(tables_dgcc) <- paste0("dgcc_ref", ref_classes)
  
  tables_drsc <- purrr::map(ref_classes, function(rc) {
    run_multinom_MI(
      data         = coverage_demographics_class,
      prob_cols    = prob_cols_drsc,
      outcome_name = "DRSC",
      n_imp        = n_imp,
      ref_class    = rc,
      zone_ref     = zone_ref,
      urban_ref    = urban_ref
    )
  })
  names(tables_drsc) <- paste0("drsc_ref", ref_classes)
  
  c(tables_dgcc, tables_drsc)
}