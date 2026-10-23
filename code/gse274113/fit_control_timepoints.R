#!/usr/bin/env Rscript

source("code/functions/runtime.R")
configure_task_runtime()
# Fit independent endpoints; never derive them from interaction coefficients.
suppressPackageStartupMessages(library(data.table))
repo <- normalizePath(getwd(), mustWork = TRUE)
source("code/functions/network_fit.R")
candidate_root <- file.path(repo, "data/gse274113/objects/GSE274113_control_full/candidates")
options <- list(input = file.path(candidate_root, "GSE274113_NT_timepoint_candidates.rds"),
                validation = file.path(candidate_root, "candidate_validation.tsv"),
                output = file.path(repo, "runs/gse274113/method_fits/multiCSN"),
                cores = 1L)
if (options$cores != 1L) stop("Only one inference core is supported.")
input <- normalizePath(options$input, mustWork = TRUE)
validation <- normalizePath(options$validation, mustWork = TRUE)
v <- fread(validation)
if (!all(c("field", "value") %in% names(v)) ||
    !identical(v[field == "training_control_cells", as.character(value)], "11447") ||
    !identical(v[field == "filter_mode", as.character(value)], "unfiltered")) {
  stop("Expected the complete, unfiltered 11,447 NT-control candidate object.")
}
out <- normalizePath(options$output, mustWork = FALSE)
if (file.exists(out)) stop("Output already exists: ", out)
stage <- paste0(out, ".staging")
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
cells <- c(day_7 = 3304L, day_9 = 4167L, day_11 = 1903L, day_14 = 2073L)
summary <- list()
for (slug in names(cells)) {
  timepoint <- sub("_", " ", slug, fixed = TRUE)
  destination <- file.path(stage, slug)
  if (file.exists(destination)) stop("Completed timepoint already exists; inspect it before resuming: ", destination)
  fit_state_networks(
    dataset = "GSE274113", input = input,
    output = destination, state = timepoint, expected_cells = cells[[slug]],
    profile = "gse274113"
  )
  endpoints <- fread(file.path(destination, "endpoint_summary.tsv"))
  summary[[slug]] <- data.table(timepoint = timepoint, training_control_cells = cells[[slug]],
    tf_gene = endpoints[endpoint == "TF-gene", selected],
    tf_region = endpoints[endpoint == "TF-region", selected],
    region_gene = endpoints[endpoint == "region-gene", selected],
    triplets = endpoints[endpoint == "factorized-chain", selected])
}
fwrite(rbindlist(summary), file.path(stage, "timepoint_network_summary.tsv"), sep = "\t")
if (!file.rename(stage, out)) stop("Atomic output promotion failed.")
cat("factorized_control_timepoints\tPASS\t", out, "\n", sep = "")
