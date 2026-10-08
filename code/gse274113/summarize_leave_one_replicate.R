#!/usr/bin/env Rscript
summarize_leave_one_replicate <- function(timepoints = "results/gse274113/replicate_deletion_stability/timepoints",
  output = "results/gse274113/replicate_deletion_stability/stability_summary") {

suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)

# ---- Parameters ----
# timepoints and output directory configuration.
values <- list(timepoints = timepoints, output = output)
# ---------------------------------------------------------------
timepoint_root <- normalizePath(values$timepoints, mustWork = TRUE)
output_root <- normalizePath(values$output, mustWork = FALSE)
if (file.exists(output_root)) stop("Output already exists: ", output_root)
slugs <- c("day_7", "day_9", "day_11", "day_14")
expected_folds <- c(day_7 = 4L, day_9 = 4L, day_11 = 3L, day_14 = 3L)
summaries <- rbindlist(lapply(slugs, function(slug) {
  path <- file.path(timepoint_root, slug, "fold_summary.tsv")
  x <- fread(path)
  if (nrow(x) != expected_folds[[slug]]) {
    stop("Unexpected fold count at ", slug)
  }
  x
}))
if (nrow(summaries) != 14L || anyDuplicated(summaries[, .(
  timepoint, heldout_replicate
)])) {
  stop("The 14-fold biological-replicate contract failed.")
}

metrics <- c(
  "edge_jaccard", "shared_edge_sign_agreement",
  "shared_edge_weight_spearman", "median_target_jaccard",
  "mean_target_jaccard", "median_regulator_jaccard",
  "mean_regulator_jaccard"
)
if (any(!vapply(summaries[, ..metrics], function(x) {
  all(is.finite(x)) && all(x >= 0 & x <= 1)
}, logical(1)))) {
  stop("A stability metric is non-finite or outside [0,1].")
}

summarize_group <- function(x, level) {
  result <- melt(
    x,
    measure.vars = metrics, variable.name = "metric", value.name = "value"
  )[, .(
    folds = .N, minimum = min(value), first_quartile = quantile(value, 0.25),
    median = median(value), third_quartile = quantile(value, 0.75),
    maximum = max(value), mean = mean(value), standard_deviation = sd(value)
  ), by = metric]
  result[, summary_level := level]
  setcolorder(result, c("summary_level", setdiff(names(result), "summary_level")))
  result
}
by_timepoint <- summaries[, summarize_group(.SD, unique(timepoint)), by = timepoint]
overall <- summarize_group(summaries, "all biological replicates")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
fwrite(summaries, file.path(output_root, "all_fold_metrics.tsv"), sep = "\t")
fwrite(by_timepoint, file.path(output_root, "stability_by_timepoint.tsv"), sep = "\t")
fwrite(overall, file.path(output_root, "stability_overall.tsv"), sep = "\t")
cat("gse274113_leave_one_replicate_summary\tPASS\t", output_root, "\n", sep = "")
}



if (sys.nframe() == 0L) summarize_leave_one_replicate()
