# scripts/fit_models.R
#
# Standalone script for fitting all hlme trajectory models.
# Run this script MANUALLY — it takes approximately 4 days.
#
# Prerequisites:
#   - targets pipeline must be up to date (run tar_make() first)
#   - Output: all_named_models.rds in the project root
#     (targets watches this file via the all_named_models_file target)
#
# Usage:
#   source("scripts/fit_models.R")   # or run interactively

library(targets)
library(lcmm)
library(future)
library(future.apply)
library(dplyr)

# -------------------------------------------------------------------------
# 1. Load coverage_data from the targets pipeline
#    This ensures models are always trained on pipeline-managed data.
# -------------------------------------------------------------------------

targets::tar_load(coverage_data)

# -------------------------------------------------------------------------
# 2. Set up parallel backend
# -------------------------------------------------------------------------

plan(multisession, workers = availableCores() - 8)
set.seed(1234)

if (!dir.exists("models")) dir.create("models")

session_start_time <- Sys.time()

# -------------------------------------------------------------------------
# 3. Define model structures
# -------------------------------------------------------------------------

model_structures <- list(
  linear_nre_homocedastic       = list(fixed = "1 + year",                          random = "~ -1",       nwg = FALSE, idiag = FALSE),
  linear_nre_heterocedastic     = list(fixed = "1 + year",                          random = "~ -1",       nwg = TRUE,  idiag = FALSE),
  quadratic_nre                 = list(fixed = "1 + year + I(year^2)",               random = "~ -1",       nwg = FALSE, idiag = FALSE),
  cubic_nre                     = list(fixed = "1 + year + I(year^2) + I(year^3)",   random = "~ -1",       nwg = FALSE, idiag = FALSE),
  linear_random_intercept       = list(fixed = "1 + year",                          random = "~ 1",        nwg = FALSE, idiag = FALSE),
  linear_random_intercept_slope = list(fixed = "1 + year",                          random = "~ 1 + year", nwg = FALSE, idiag = FALSE),
  quadratic_random_effects      = list(fixed = "1 + year + I(year^2)",               random = "~ 1 + year", nwg = FALSE, idiag = FALSE),
  quadratic_random_effects_prop = list(fixed = "1 + year + I(year^2)",               random = "~ 1 + year", nwg = TRUE,  idiag = FALSE),
  cubic_random_effects          = list(fixed = "1 + year + I(year^2) + I(year^3)",   random = "~ 1 + year", nwg = FALSE, idiag = FALSE),
  cubic_random_effects_prop     = list(fixed = "1 + year + I(year^2) + I(year^3)",   random = "~ 1 + year", nwg = TRUE,  idiag = FALSE)
)

# -------------------------------------------------------------------------
# 4. Build task list: all combinations of structure × outcome × n_classes
# -------------------------------------------------------------------------

model_tasks <- list()

for (structure_name in names(model_structures)) {
  structure <- model_structures[[structure_name]]
  for (var in c("drsc", "dgcc")) {
    for (ng in 1:7) {
      key <- paste(structure_name, var, ng, sep = ":")
      model_tasks[[key]] <- list(
        var            = var,
        structure      = structure,
        ng             = ng,
        structure_name = structure_name
      )
    }
  }
}

# -------------------------------------------------------------------------
# 5. Fit models in parallel
#    Each model is saved individually to models/ as it completes.
#    If a model fails, the error is logged and the loop continues.
# -------------------------------------------------------------------------

results <- future_lapply(model_tasks, function(task) {

  set.seed(1234)

  var            <- task$var
  structure      <- task$structure
  ng             <- task$ng
  structure_name <- task$structure_name
  model_name     <- paste0(ng, "class_", structure_name, "_", var, "_model")
  rds_path       <- file.path("models", paste0(model_name, ".rds"))

  # Skip if already fitted (useful for resuming after interruption)
  if (file.exists(rds_path)) {
    message("⏭️  Skipping (already exists): ", model_name)
    return(model_name)
  }

  fixed_formula   <- as.formula(paste0(var, " ~ ", structure$fixed))
  random_formula  <- as.formula(structure$random)
  mixture_formula <- if (ng > 1) fixed_formula else NULL
  nwg_value       <- if (ng > 1) structure$nwg else FALSE

  tryCatch({
    if (ng == 1) {
      fitted_model <- lcmm::hlme(
        fixed   = fixed_formula,
        random  = random_formula,
        subject = "id",
        ng      = ng,
        nwg     = nwg_value,
        idiag   = structure$idiag,
        data    = coverage_data
      )
    } else {
      one_class_model <- lcmm::hlme(
        fixed   = fixed_formula,
        random  = random_formula,
        subject = "id",
        ng      = 1,
        nwg     = FALSE,
        idiag   = structure$idiag,
        data    = coverage_data
      )
      fitted_model <- lcmm::gridsearch(
        rep     = 20,
        maxiter = 1000,
        minit   = one_class_model,
        lcmm::hlme(
          fixed   = fixed_formula,
          mixture = mixture_formula,
          random  = random_formula,
          subject = "id",
          ng      = ng,
          nwg     = nwg_value,
          idiag   = structure$idiag,
          data    = coverage_data
        )
      )
    }

    saveRDS(fitted_model, file = rds_path)
    message("✅ Success: ", model_name)
    return(model_name)

  }, error = function(e) {
    message("❌ Failed: ", model_name, " — ", e$message)
    return(NULL)
  })
})

message("Total elapsed: ", round(difftime(Sys.time(), session_start_time, units = "hours"), 1), " hours")

# -------------------------------------------------------------------------
# 6. Collect all saved models into a single named list and save
#    Saving to the project root so targets can track it via all_named_models_file.
# -------------------------------------------------------------------------

model_files <- list.files("models", pattern = "\\.rds$", full.names = TRUE)

all_named_models <- list()
for (file in model_files) {
  name <- tools::file_path_sans_ext(basename(file))
  all_named_models[[name]] <- readRDS(file)
}

saveRDS(all_named_models, "all_named_models.rds")

message("Done. ", length(all_named_models), " models saved to all_named_models.rds")

# -------------------------------------------------------------------------
# After this script completes, run tar_make() to propagate the updated
# all_named_models.rds through the rest of the targets pipeline.
# -------------------------------------------------------------------------
