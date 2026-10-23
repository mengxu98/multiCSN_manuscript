#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
source("code/functions/network_io.R")
root <- "results/gse274113/networks"
methods <- c("multiCSN", "Pando", "SCENICplus")
if (!all(methods %in% list.dirs(root, recursive = FALSE, full.names = FALSE)))
  stop("Missing canonical method suite.")
endpoints <- c("tf_gene", "tf_region", "region_gene", "triplet")
rows <- list()
for (method in methods) for (slug in c("day_7", "day_9", "day_11", "day_14")) {
  for (endpoint in endpoints) {
    path <- file.path(root, method, slug, paste0(endpoint, ".tsv.gz"))
    if (!file.exists(path)) stop("Missing canonical endpoint: ", path)
    rows[[length(rows) + 1L]] <- data.table(method = method, timepoint = slug,
      endpoint = endpoint, relations = nrow(read_fit_checkpoint(path)),
      source = "canonical endpoint table")
  }
}
summary <- rbindlist(rows)
dir.create(root, recursive = TRUE, showWarnings = FALSE)
write_tsv(summary, file.path(root, "endpoint_summary.tsv"))
cat("Catalog derived from", nrow(summary), "endpoint tables\n")
