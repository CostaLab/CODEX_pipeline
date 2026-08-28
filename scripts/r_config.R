# Shared config loader for all R steps of the CODEX pipeline.
#
# Usage in a script:
#   args <- commandArgs(trailingOnly = TRUE)
#   config_path <- if (length(args) >= 1) args[1] else file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), "..", "config.yaml")
#   source(file.path(dirname(config_path_of_this_file), "r_config.R"))
#   cfg <- load_config(config_path)
#   pth <- get_pipeline_paths(cfg)

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required (install.packages('yaml')).")
}
library(yaml)

load_config <- function(config_path) {
  yaml::read_yaml(config_path)
}

get_pipeline_paths <- function(cfg) {
  wdir <- cfg$paths$wdir
  pipeline_run <- file.path(wdir, cfg$paths$pipeline_run_subdir)
  list(
    wdir = wdir,
    pipeline_run = pipeline_run,
    quant_dir = file.path(pipeline_run, cfg$paths$quant_subdir),
    visual_output_dir = file.path(pipeline_run, cfg$paths$visual_output_subdir),
    out_data_dir = paste0(pipeline_run, "/"),
    fig_dir = file.path(pipeline_run, cfg$paths$visual_output_subdir)
  )
}

# Build named lists keyed by "run<id>", matching what the original scripts
# hardcoded (condition_vecs, working_samples).
get_run_lists <- function(cfg) {
  condition_vecs <- list()
  working_samples <- list()
  for (run in cfg$runs) {
    key <- paste0("run", run$id)
    condition_vecs[[key]] <- unlist(run$conditions)
    ws <- run$working_samples
    if (is.null(ws)) ws <- seq_along(run$conditions)
    working_samples[[key]] <- unlist(ws)
  }
  list(condition_vecs = condition_vecs, working_samples = working_samples)
}

get_config_path_from_args <- function(default_relative = "../config.yaml") {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) {
    return(args[1])
  }
  # fall back to config.yaml next to this script's parent dir
  this_file <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(this_file) == 0) return(default_relative)
  file.path(dirname(this_file), default_relative)
}
