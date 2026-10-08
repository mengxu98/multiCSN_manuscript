#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
source("code/functions/network_io.R")
args <- list(input = "results/gse274113/networks",
                            output = "results/gse274113/networks")
# LINGER currently has only retained run summaries, not these canonical tables.
# Never infer its endpoint availability from directory existence or summary counts.
methods <- c("multiCSN", "Pando", "SCENICplus")
if (!all(methods %in% list.dirs(args$input, recursive = FALSE, full.names = FALSE)))
  stop("Missing canonical method suite.")
endpoints <- c("tf_gene", "tf_region", "region_gene", "triplet")
rows <- list()
for (method in methods) for (slug in c("day_7", "day_9", "day_11", "day_14")) {
  for (endpoint in endpoints) {
    path <- file.path(args$input, method, slug, paste0(endpoint, ".tsv.gz"))
    if (!file.exists(path)) stop("Missing canonical endpoint: ", path)
    rows[[length(rows) + 1L]] <- data.table(method = method, timepoint = slug,
      endpoint = endpoint, relations = nrow(read_fit_checkpoint(path)),
      source = "canonical endpoint table; historical provenance is not inferred")
  }
}
summary <- rbindlist(rows)
availability <- unique(summary[, .(method, endpoint)])
availability[, `:=`(available = TRUE, reason = "endpoint tables present at all four timepoints")]
paths <- file.path(args$output, c("endpoint_summary.tsv", "endpoint_availability.tsv"))
if (any(file.exists(paths))) stop("Archive existing catalog tables before replacement.")
dir.create(args$output, recursive = TRUE, showWarnings = FALSE)
write_tsv(summary, paths[[1]])
write_tsv(availability, paths[[2]])
cat("PASS catalog derived from", nrow(summary), "canonical endpoint tables\n")
