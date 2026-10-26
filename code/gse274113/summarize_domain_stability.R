#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)

# ---- Parameters ----
# run_root: leave-one-replicate tree; output: summary output directory.
arguments <- list(
  run_root = "results/gse274113/stability",
  output = "results/gse274113/stability"
)
# ---------------------------------------------------------------
rebuilt_root <- normalizePath(arguments$run_root, mustWork = TRUE)
input_root <- rebuilt_root
output_root <- normalizePath(arguments$output, mustWork = FALSE)

input_paths <- c(
  folds = file.path(input_root, "fold_decomposition.tsv"),
  pairwise = file.path(
    rebuilt_root,
    "pairwise_fold_metrics.tsv"
  ),
  top_rank = file.path(input_root, "top_rank_decomposition.tsv"),
  rank_decile = file.path(input_root, "rank_decile_decomposition.tsv"),
  layer_attrition = file.path(input_root, "layer_attrition.tsv")
)

folds <- fread(input_paths[["folds"]])
pairwise <- fread(input_paths[["pairwise"]])
top_rank <- fread(input_paths[["top_rank"]])
rank_decile <- fread(input_paths[["rank_decile"]])
layer_attrition <- fread(input_paths[["layer_attrition"]])

timepoint_order <- c("day 7", "day 9", "day 11", "day 14")
observed_folds <- folds[, .(N = .N,
  complete_replicate_count = unique(training_replicates + 1L)), by = timepoint]
if (anyNA(folds[, .(timepoint, heldout_replicate, training_replicates)]) ||
  anyDuplicated(folds[, .(timepoint, heldout_replicate)]) ||
  !setequal(observed_folds$timepoint, timepoint_order) || any(observed_folds$N < 2L) ||
  anyDuplicated(observed_folds$timepoint) ||
  any(observed_folds$N != observed_folds$complete_replicate_count)) {
  stop("Invalid replicate metadata.")
}
expected_pairs <- setNames(choose(observed_folds$N, 2L), observed_folds$timepoint)[timepoint_order]
observed_pairs <- pairwise[, .N, by = timepoint]
pairwise_metrics <- c(
  "edge_jaccard", "shared_edge_sign_agreement",
  "shared_edge_weight_spearman"
)
if (nrow(pairwise) != sum(expected_pairs) ||
  anyDuplicated(pairwise[, .(timepoint, fold_left, fold_right)]) ||
  !setequal(observed_pairs$timepoint, timepoint_order) ||
  any(setNames(observed_pairs$N, observed_pairs$timepoint)[timepoint_order] != expected_pairs) ||
  !all(pairwise_metrics %in% names(pairwise)) ||
  any(!vapply(pairwise[, ..pairwise_metrics], function(x) {
    all(is.finite(x)) && all(x >= 0 & x <= 1)
  }, logical(1)))) {
  stop("The deletion-fold pair design failed.")
}

required_fraction_columns <- c(
  "rna_model_fraction", "candidate_retention_given_rna_model",
  "fit_eligibility_given_rna_and_candidate", "selection_given_fit_eligible",
  "unconditional_selection_fraction"
)
if (!all(required_fraction_columns %in% names(top_rank)) ||
  !all(required_fraction_columns %in% names(rank_decile))) {
  stop("Decomposition fractions are incomplete.")
}
for (column in required_fraction_columns) {
  if (any(top_rank[[column]] < 0 | top_rank[[column]] > 1, na.rm = TRUE) ||
    any(rank_decile[[column]] < 0 | rank_decile[[column]] > 1,
      na.rm = TRUE
    )) {
    stop("Out-of-range fraction: ", column)
  }
}
fraction_denominators <- c(
  rna_model_fraction = "fold_opportunities",
  candidate_retention_given_rna_model = "rna_model_eligible",
  fit_eligibility_given_rna_and_candidate = "candidate_and_rna_model",
  selection_given_fit_eligible = "fit_eligible",
  unconditional_selection_fraction = "fold_opportunities"
)
for (column in names(fraction_denominators)) {
  denominator <- fraction_denominators[[column]]
  if (any(is.na(top_rank[[column]]) & top_rank[[denominator]] > 0) ||
    any(is.na(rank_decile[[column]]) & rank_decile[[denominator]] > 0)) {
    stop("Missing fraction with a non-zero denominator: ", column)
  }
}

quantile_value <- function(x, probability) {
  as.numeric(stats::quantile(x, probability,
    na.rm = TRUE, names = FALSE,
    type = 7L
  ))
}

timepoint_summary <- folds[, .(
  folds = .N,
  full_edges = unique(full_edges),
  training_cells_minimum = min(training_cells),
  training_cells_maximum = max(training_cells),
  fold_edges_minimum = min(fold_edges),
  median_fold_edges = as.numeric(stats::median(fold_edges)),
  fold_edges_maximum = max(fold_edges),
  elapsed_seconds_minimum = min(elapsed_seconds),
  elapsed_seconds_median = as.numeric(stats::median(elapsed_seconds)),
  elapsed_seconds_maximum = max(elapsed_seconds),
  edge_jaccard_q1 = quantile_value(edge_jaccard, 0.25),
  edge_jaccard_median = as.numeric(stats::median(edge_jaccard)),
  edge_jaccard_q3 = quantile_value(edge_jaccard, 0.75),
  sign_agreement_minimum = min(shared_edge_sign_agreement),
  sign_agreement_median = as.numeric(stats::median(shared_edge_sign_agreement)),
  weight_spearman_q1 = quantile_value(shared_edge_weight_spearman, 0.25),
  weight_spearman_median = as.numeric(stats::median(shared_edge_weight_spearman)),
  weight_spearman_q3 = quantile_value(shared_edge_weight_spearman, 0.75),
  rna_model_fraction_median = as.numeric(stats::median(
    full_edge_rna_model_eligible / full_edges
  )),
  candidate_retention_given_rna_model_median = as.numeric(stats::median(
    full_edge_candidate_retention_given_rna_model
  )),
  fit_eligibility_given_rna_and_candidate_median = as.numeric(stats::median(
    full_edge_fit_eligibility_given_rna_and_candidate
  )),
  selection_given_fit_eligible_median = as.numeric(stats::median(
    full_edge_selection_given_fit_eligible
  )),
  unconditional_selection_fraction_median = as.numeric(stats::median(
    full_edge_selected / full_edges
  ))
), by = timepoint]
timepoint_summary[, timepoint_index := match(timepoint, timepoint_order)]
setorder(timepoint_summary, timepoint_index)
timepoint_summary[, timepoint_index := NULL]

pairwise_summary <- pairwise[, .(
  fold_pairs = .N,
  edge_jaccard_q1 = quantile_value(edge_jaccard, 0.25),
  edge_jaccard_median = as.numeric(stats::median(edge_jaccard)),
  edge_jaccard_q3 = quantile_value(edge_jaccard, 0.75),
  sign_agreement_minimum = min(shared_edge_sign_agreement),
  sign_agreement_median = as.numeric(stats::median(shared_edge_sign_agreement)),
  weight_spearman_q1 = quantile_value(shared_edge_weight_spearman, 0.25),
  weight_spearman_median = as.numeric(stats::median(shared_edge_weight_spearman)),
  weight_spearman_q3 = quantile_value(shared_edge_weight_spearman, 0.75)
), by = timepoint]
pairwise_summary[, timepoint_index := match(timepoint, timepoint_order)]
setorder(pairwise_summary, timepoint_index)
pairwise_summary[, timepoint_index := NULL]

top_limits <- top_rank[stratum == "top_k", .(
  maximum_available = max(as.integer(stratum_value))
), by = timepoint]
top_limit <- min(100L, min(top_limits$maximum_available))
top100 <- top_rank[
  stratum == "top_k" & as.integer(stratum_value) == top_limit
]
if (nrow(top100) != length(timepoint_order) || !setequal(top100$timepoint, timepoint_order)) {
  stop("Missing common top-rank decomposition.")
}
top100[, timepoint_index := match(timepoint, timepoint_order)]
setorder(top100, timepoint_index)
top100[, timepoint_index := NULL]
rank_decile[, `:=`(
  timepoint_index = match(timepoint, timepoint_order),
  decile_index = as.integer(stratum_value)
)]
setorder(rank_decile, timepoint_index, decile_index)
rank_decile[, c("timepoint_index", "decile_index") := NULL]
decile_trends <- rank_decile[, .(
  conditional_selection_decile_spearman = suppressWarnings(stats::cor(
    as.integer(stratum_value), selection_given_fit_eligible,
    method = "spearman", use = "complete.obs"
  )),
  unconditional_selection_decile_spearman = suppressWarnings(stats::cor(
    as.integer(stratum_value), unconditional_selection_fraction,
    method = "spearman", use = "complete.obs"
  ))
), by = timepoint]
if (nrow(decile_trends) != length(timepoint_order)) {
  stop("Rank-decile stability trends are incomplete.")
}
timepoint_summary <- merge(timepoint_summary, decile_trends,
  by = "timepoint",
  sort = FALSE
)
timepoint_summary[, timepoint_index := match(timepoint, timepoint_order)]
setorder(timepoint_summary, timepoint_index)
timepoint_summary[, timepoint_index := NULL]
layer_attrition[, timepoint_index := match(timepoint, timepoint_order)]
setorder(layer_attrition, timepoint_index, heldout_replicate)
layer_attrition[, timepoint_index := NULL]

stage_root <- tempfile(".result-files-", tmpdir = dirname(output_root))
dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
dir.create(stage_root, recursive = FALSE)

fwrite(timepoint_summary, file.path(stage_root, "timepoint_stability_summary.tsv"), sep = "\t")
fwrite(pairwise_summary, file.path(stage_root, "pairwise_stability_summary.tsv"), sep = "\t")
fwrite(top100, file.path(stage_root, "top100_decomposition.tsv"), sep = "\t")

publish_result_files(stage_root, output_root)
