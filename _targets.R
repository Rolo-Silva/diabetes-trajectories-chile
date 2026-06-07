# _targets.R
# Pipeline for Chilean municipal diabetes-care trajectories analysis (2011-2023)

library(targets)
library(tarchetypes)

tar_option_set(
  packages = c(
    "targets", "dplyr", "tidyr", "readr", "stringr", "stringi",
    "lubridate", "forcats",
    "ggplot2", "glue", "purrr", "janitor", "readxl",
    "lcmm", "LCTMtools",
    "patchwork", "cowplot",
    "ggsci", "scales", "psych", "grid", "MASS",
    "tibble", "rlang",
    "nnet", "broom",
    "lme4", "lmerTest"
  )
)

# Source all function files in R/
tar_source("R")

base_dir <- "/Users/rolo/Documents/trajectories"

list(
  
  # -------------------------------------------------------------------------
  # 1. Load raw annual series (2011-2023)
  # -------------------------------------------------------------------------
  
  tar_target(SerieP2011,
             read.csv(file.path(base_dir, "SerieP2011.txt"), sep = ";") |>
               janitor::clean_names() |>
               dplyr::mutate(
                 id_establecimiento = as.numeric(stringr::str_remove_all(id_establecimiento, "-")),
                 id_establecimiento = sprintf("1%05d", id_establecimiento)
               ) |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  tar_target(SerieP2012,
             read.csv(file.path(base_dir, "SerieP2012.txt"), sep = ";") |>
               janitor::clean_names() |>
               dplyr::mutate(
                 id_establecimiento = as.numeric(stringr::str_remove_all(id_establecimiento, "-")),
                 id_establecimiento = sprintf("1%05d", id_establecimiento)
               ) |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  tar_target(SerieP2013,
             read.csv(file.path(base_dir, "SerieP2013.txt"), sep = ";") |>
               janitor::clean_names() |>
               dplyr::mutate(
                 id_establecimiento = as.numeric(stringr::str_remove_all(id_establecimiento, "-")),
                 id_establecimiento = sprintf("1%05d", id_establecimiento)
               ) |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  tar_target(SerieP2014,
             read.csv(file.path(base_dir, "SerieP2014.txt"), sep = ";") |>
               janitor::clean_names() |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  tar_target(SerieP2015,
             read.csv(file.path(base_dir, "SerieP2015.txt"), sep = ";") |>
               janitor::clean_names() |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  tar_target(SerieP2016,
             read.csv(file.path(base_dir, "SerieP2016.txt"), sep = ";") |>
               janitor::clean_names() |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  tar_target(SerieP2017_pre,
             read.csv(file.path(base_dir, "SerieP2017.txt"), sep = ";") |>
               janitor::clean_names() |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  tar_target(SerieP2018,
             read.csv(file.path(base_dir, "SerieP2018.txt"), sep = ";") |>
               janitor::clean_names() |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  tar_target(SerieP2019,
             read.csv(file.path(base_dir, "SerieP2019.txt"), sep = ";") |>
               janitor::clean_names() |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  tar_target(SerieP2020,
             read.csv(file.path(base_dir, "SerieP2020.txt"), sep = ";") |>
               janitor::clean_names() |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  tar_target(SerieP2021,
             read.csv(file.path(base_dir, "SerieP2021.txt"), sep = ";") |>
               janitor::clean_names() |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  tar_target(SerieP2022,
             read.csv(file.path(base_dir, "SerieP2022.txt"), sep = ";") |>
               janitor::clean_names() |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  tar_target(SerieP2023,
             read.csv(file.path(base_dir, "SerieP2023.txt"), sep = ";") |>
               janitor::clean_names() |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  # -------------------------------------------------------------------------
  # 2. Patch for SerieP2017
  # -------------------------------------------------------------------------
  
  tar_target(datosf_2017,
             readxl::read_excel(file.path(base_dir, "datosf_2017.xlsx")) |>
               tidyr::separate(
                 col = id_establecimiento,
                 into = c("id_establecimiento", "consultorio"),
                 sep = "-"
               ) |>
               dplyr::select(-consultorio) |>
               tidyr::gather(
                 -id_establecimiento, -id_servicio, -id_comuna, -id_region,
                 key = codigo_prestacion,
                 value = col01,
                 2:3
               ) |>
               dplyr::mutate(
                 mes = rep(12, 282),
                 ano = rep(2017, 282)
               ) |>
               dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  ),
  
  tar_target(SerieP2017,
             dplyr::bind_rows(SerieP2017_pre, datosf_2017)
  ),
  
  # -------------------------------------------------------------------------
  # 3. DPA2018 lookup table
  # -------------------------------------------------------------------------
  
  tar_target(dpa2018,
             readxl::read_excel(file.path(base_dir, "DPA2018.xls"), skip = 1) |>
               dplyr::rename(
                 comuna         = `Nombre Comuna`,
                 id_comuna      = `Código Comuna desde 2010`,
                 id_comuna2     = `Código Comuna desde 2018`,
                 region         = `Nombre Región desde 2008`,
                 region2        = `Nombre Región desde 2018`,
                 id_region      = `Código Región desde 2008`,
                 id_region2     = `Código Región desde 2018`,
                 servicio_salud = `Nombre Servicio de Salud desde 2008`,
                 id_servicio2   = `Código Servicio Salud desde 2008`
               ) |>
               dplyr::mutate(
                 servicio_salud = dplyr::case_when(
                   servicio_salud == "Viña Del Mar Quillota" ~ "Viña del Mar Quillota",
                   TRUE ~ servicio_salud
                 ),
                 id_region = as.numeric(id_region)
               ) |>
               dplyr::slice(-dplyr::n())
  ),
  
  # -------------------------------------------------------------------------
  # 4. Combine all series → filter for diabetes → build coverage_data
  # -------------------------------------------------------------------------
  
  tar_target(all_series,
             make_all_series(list(
               SerieP2011, SerieP2012, SerieP2013, SerieP2014, SerieP2015, SerieP2016,
               SerieP2017, SerieP2018, SerieP2019, SerieP2020, SerieP2021, SerieP2022,
               SerieP2023
             ))
  ),
  
  tar_target(all_series_diabetes,
             make_all_series_diabetes(all_series, dpa2018)
  ),
  
  tar_target(coverage_data,
             make_coverage_data(all_series_diabetes)
  ),
  
  tar_target(idc,
             readxl::read_excel(file.path(base_dir, "indice_desarrollo_comunal_pdf_excel.xlsx")) |>
               janitor::clean_names()
  ),
  
  tar_target(comunas_pndr,
             readxl::read_excel(file.path(base_dir, "comunas_pndr.xlsx")) |>
               janitor::clean_names()
  ),
  
  # -------------------------------------------------------------------------
  # 5. Models — fitted externally (see scripts/fit_models.R)
  # -------------------------------------------------------------------------
  
  tar_target(
    all_named_models_file,
    "all_named_models.rds",
    format = "file"
  ),
  
  tar_target(
    all_named_models,
    {
      obj <- readRDS(all_named_models_file)
      if (!is.list(obj))
        stop("all_named_models.rds must contain a list of fitted lcmm models.")
      if (is.null(names(obj)))
        stop("all_named_models.rds must be a named list. names(obj) is NULL.")
      
      required_models <- c(
        "4class_cubic_nre_dgcc_model",
        "4class_cubic_nre_drsc_model",
        "4class_linear_nre_homocedastic_dgcc_model",
        "4class_linear_nre_homocedastic_drsc_model"
      )
      missing_models <- setdiff(required_models, names(obj))
      if (length(missing_models) > 0)
        stop(paste0("Missing required models in all_named_models.rds: ",
                    paste(missing_models, collapse = ", ")))
      obj
    }
  ),
  
  tar_target(
    model_inventory,
    tibble::tibble(
      model_name  = names(all_named_models),
      model_class = purrr::map_chr(all_named_models, ~ paste(class(.x), collapse = "; ")),
      converged   = purrr::map_int(all_named_models,
                                   ~ if (is.null(.x$conv)) NA_integer_ else .x$conv)
    )
  ),
  
  # -------------------------------------------------------------------------
  # 6. Model adequacy table
  # -------------------------------------------------------------------------
  
  tar_target(
    model_adequacy_table,
    make_model_adequacy_table(all_named_models, coverage_data)
  ),
  
  # -------------------------------------------------------------------------
  # 7. BIC elbow plots — Figure 2
  # -------------------------------------------------------------------------
  
  tar_target(
    bic_elbow_plot,
    make_bic_elbow_plot(
      model_adequacy_table,
      out_path = "outputs/figures/fig2_bic_elbow.png",
      width = 12, height = 8, dpi = 300
    ),
    format = "file"
  ),
  
  # -------------------------------------------------------------------------
  # 8. Predicted mean trajectory plots — Figure 3
  # -------------------------------------------------------------------------
  
  tar_target(
    predicted_means_plot,
    make_predicted_means_plot(
      all_named_models,
      out_path = "outputs/figures/fig3_predicted_means.png",
      width = 12, height = 8, dpi = 300
    ),
    format = "file"
  ),
  
  # -------------------------------------------------------------------------
  # 9. Spaghetti plots — Figure 4
  # -------------------------------------------------------------------------
  
  tar_target(
    spaghetti_plot,
    make_spaghetti_plot(
      all_named_models,
      trajectory_labels = list(
        "4class_cubic_nre_dgcc_model" = c("Lowest", "Highest", "Stable upper medium", "Stable medium"),
        "4class_cubic_nre_drsc_model" = c("Increasing", "Stable medium", "Predominantly Highest", "Lowest")
      ),
      out_path = "outputs/figures/fig4_spaghetti.png",
      width = 12, height = 7, dpi = 300
    ),
    format = "file"
  ),
  
  # -------------------------------------------------------------------------
  # 10. Elsensohn residual envelope plots — Figure 5
  # -------------------------------------------------------------------------
  
  tar_target(
    elsensohn_plot,
    make_elsensohn_plot(
      all_named_models,
      trajectory_labels = list(
        "4class_cubic_nre_dgcc_model" = c("Lowest", "Highest", "Stable upper medium", "Stable medium"),
        "4class_cubic_nre_drsc_model" = c("Increasing", "Stable medium", "Predominantly Highest", "Lowest")
      ),
      out_path = "outputs/figures/fig5_elsensohn.png",
      width = 12, height = 5, dpi = 300
    ),
    format = "file"
  ),
  
  # -------------------------------------------------------------------------
  # 11. Standardised residual benchmark plots — Figure 6
  # -------------------------------------------------------------------------
  
  tar_target(
    residual_benchmark_plot,
    make_residual_benchmark_plot(
      all_named_models,
      out_path = "outputs/figures/fig6_residual_benchmark.png",
      width = 12, height = 7, dpi = 300
    ),
    format = "file"
  ),
  
  # -------------------------------------------------------------------------
  # 12. Enrichment: coverage_data + IDC + PNDR + classes
  # -------------------------------------------------------------------------
  
  tar_target(
    coverage_data_enriched,
    make_coverage_data_enriched(coverage_data, idc, comunas_pndr)
  ),
  
  tar_target(
    coverage_classes,
    make_coverage_classes(all_named_models)
  ),
  
  tar_target(
    coverage_demographics_class,
    make_coverage_demographics_class(coverage_data_enriched, coverage_classes)
  ),
  
  # -------------------------------------------------------------------------
  # 13. Trajectory plots for manuscript — SF3 + Figure 1B/2B
  # -------------------------------------------------------------------------
  
  tar_target(
    sf3_mean_trends,
    make_mean_trends_plot(
      coverage_data = coverage_data,
      out_path      = "outputs/figures/sf3_mean_trends.png"
    ),
    format = "file"
  ),
  
  tar_target(
    coverage_data_with_classes,
    coverage_data |> dplyr::left_join(coverage_classes, by = "id")
  ),
  
  tar_target(
    figure_1b_egcc_trajectories,
    make_class_trajectories_plot(
      coverage_data_with_classes = coverage_data_with_classes,
      outcome      = dgcc,
      class_var    = class_dgcc,
      class_labels = c(
        "Lowest: n=26 (8.7%)",
        "Highest: n=27 (9.0%)",
        "Stable upper-medium: n=134 (44.8%)",
        "Stable medium: n=112 (37.5%)"
      ),
      y_axis_label = "Municipality mean EGCC",
      out_path     = "outputs/figures/figure_1b_egcc_trajectories.png"
    ),
    format = "file"
  ),
  
  tar_target(
    figure_2b_drsc_trajectories,
    make_class_trajectories_plot(
      coverage_data_with_classes = coverage_data_with_classes,
      outcome      = drsc,
      class_var    = class_drsc,
      class_labels = c(
        "Increasing: n=43 (14.4%)",
        "Stable medium: n=97 (32.4%)",
        "Predominantly highest: n=58 (19.4%)",
        "Lowest: n=101 (33.8%)"
      ),
      y_axis_label = "Municipality mean DRSC",
      out_path     = "outputs/figures/figure_2b_drsc_trajectories.png"
    ),
    format = "file"
  ),
  
  # -------------------------------------------------------------------------
  # 14. Sensitivity validation: Stage 1 + Stage 2 (Figures 7 + 8)
  # Default: 1 seed each (fast pipeline).
  # For production: change `sensitivity_seeds` to 1:20 (or any vector).
  # -------------------------------------------------------------------------
  
  tar_target(sensitivity_seeds, 1L),
  
  tar_target(
    subsampling_validation,
    make_subsampling_validation(
      coverage_data    = coverage_data,
      all_named_models = all_named_models,
      seed_vec         = sensitivity_seeds,
      out_path         = "outputs/figures/fig7_subsampling.png"
    )
  ),
  
  tar_target(
    perturbation_validation,
    make_perturbation_validation(
      coverage_data    = coverage_data,
      all_named_models = all_named_models,
      seed_vec         = sensitivity_seeds,
      out_path         = "outputs/figures/fig8_perturbation.png"
    )
  ),
  
  tar_target(
    multinom_validation,
    make_multinom_validation(
      coverage_demographics_class = coverage_demographics_class,
      zone_ref  = "Centre",
      urban_ref = "Mixed",
      n_imp     = 10
    )
  ),
  
  # -------------------------------------------------------------------------
  # 16. Cross-classification: Supplementary Figure 9
  # 4×4 contingency table of EGCC × DRSC trajectory classes with
  # chi-square (Monte Carlo), Cramér's V, and heatmap figure.
  # -------------------------------------------------------------------------
  
  tar_target(
    cross_classification,
    make_cross_classification(
      coverage_demographics_class = coverage_demographics_class,
      out_path  = "outputs/figures/sf9_cross_classification.png",
      chi_sim_B = 10000
    )
  ),
  
  # -------------------------------------------------------------------------
  # 17. Mixed-effects analysis: Supplementary Table 7
  # Fits MAIN and INTERACTION lmer models for both outcomes; returns the
  # four fits, tidy coefficient tables with Wald 95% CIs, LRT comparison
  # (main vs interaction), and preformatted Supplementary Table 7.
  # -------------------------------------------------------------------------
  
  tar_target(
    mixed_effects_analysis,
    make_mixed_effects_analysis(
      coverage_data_enriched = coverage_data_enriched
    )
  )
  
)

