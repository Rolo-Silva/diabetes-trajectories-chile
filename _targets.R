# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)

tar_option_set(
  packages = c(
    "targets", "dplyr", "tidyr", "readr", "stringr", "lubridate", "forcats",
    "ggplot2", "glue", "purrr", "janitor", "readxl", "lcmm", "LCTMtoolkit"
  )
)

# Load all functions from Utils.R
tar_source("R/Utils.R")

base_dir <- "/Users/rolo/Documents/trajectories"

list(

  # -------------------------------------------------------------------------
  # 1. Load raw annual series (2011-2023)
  #    All columns cast to character for safe row-binding across years.
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
  #    Uses dplyr::bind_rows() — no plyr dependency needed.
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
        key = codigo_prestacion, value = col01, 2:3
      ) |>
      dplyr::mutate(mes = rep(12, 282), ano = rep(2017, 282)) |>
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
      dplyr::slice(-dplyr::n())   # remove trailing summary row
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

  # -------------------------------------------------------------------------
  # 5. Models — fitted externally (takes ~4 days).
  #    targets tracks the .rds file; if it changes, downstream targets rerun.
  # -------------------------------------------------------------------------

  tar_target(all_named_models_file,
    "all_named_models.rds",
    format = "file"
  ),

  tar_target(all_named_models,
    readRDS(all_named_models_file)
  ),

  # -------------------------------------------------------------------------
  # 6. Model adequacy table — built from the loaded models
  # -------------------------------------------------------------------------

  tar_target(model_adequacy_table,
    make_model_adequacy_table(all_named_models)
  )

)
