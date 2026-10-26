#!/usr/bin/env Rscript

source("code/functions/runtime.R")
configure_task_runtime()
# Fit independent endpoints; never derive them from interaction coefficients.
suppressPackageStartupMessages(library(data.table))
repo <- normalizePath(getwd(), mustWork = TRUE)
source("code/functions/network_fit.R")
candidate_root <- file.path(repo, "data/gse274113/objects/GSE274113_control_full/candidates")
options <- list(input = file.path(candidate_root, "GSE274113_NT_timepoint_candidates.rds"),
                output = file.path(repo, "runs/gse274113/method_fits/multiCSN"),
                cores = 1L)
if (options$cores != 1L) stop("Only one inference core is supported.")
input <- normalizePath(options$input, mustWork = TRUE)
out <- normalizePath(options$output, mustWork = FALSE)
if (file.exists(out)) stop("Output already exists: ", out)
stage <- paste0(out, ".staging")
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
timepoints <- c(day_7 = "day 7", day_9 = "day 9", day_11 = "day 11", day_14 = "day 14")
summary <- list()
for (slug in names(timepoints)) {
  timepoint <- timepoints[[slug]]
  destination <- file.path(stage, slug)
  if (file.exists(destination)) stop("Completed timepoint already exists; inspect it before resuming: ", destination)
  fitted_cells <- fit_state_networks(
    dataset = "GSE274113", input = input,
    output = destination, state = timepoint,
    profile = "gse274113"
  )
  endpoints <- fread(file.path(destination, "endpoint_summary.tsv"))
  summary[[slug]] <- data.table(timepoint = timepoint, training_control_cells = nrow(fitted_cells),
    tf_gene = endpoints[endpoint == "TF-gene", selected],
    tf_region = endpoints[endpoint == "TF-region", selected],
    region_gene = endpoints[endpoint == "region-gene", selected],
    triplets = endpoints[endpoint == "factorized-chain", selected])
}
fwrite(rbindlist(summary), file.path(stage, "timepoint_network_summary.tsv"), sep = "\t")
if (!file.rename(stage, out)) stop("Atomic output promotion failed.")
