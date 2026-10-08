#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)


safe_fraction <- function(numerator, denominator) {
  fifelse(denominator > 0, numerator / denominator, NA_real_)
}

# ---- Parameters ----
# timepoints: replicate_deletion_stability timepoints directory.
arguments <- list(
  timepoints = "results/gse274113/replicate_deletion_stability/timepoints",
  output = "results/gse274113/replicate_deletion_stability/multilayer_decomposition"
)
# ---------------------------------------------------------------
timepoint_root <- normalizePath(arguments$timepoints, mustWork = TRUE)
output_root <- normalizePath(arguments$output, mustWork = FALSE)
if (file.exists(output_root)) stop("Output already exists: ", output_root)
dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
stage_root <- paste0(output_root, ".staging_", format(
  Sys.time(), "%Y%m%dT%H%M%SZ",
  tz = "UTC"
))
dir.create(stage_root, recursive = FALSE)

timepoint_lookup <- c(
  day_7 = "day 7", day_9 = "day 9", day_11 = "day 11", day_14 = "day 14"
)
expected_folds <- c(day_7 = 4L, day_9 = 4L, day_11 = 3L, day_14 = 3L)
required_state_columns <- c(
  "rna_detection_eligible", "rna_model_eligible",
  "multiomic_candidate_present", "fit_eligible", "selected_in_fold"
)
fold_rows <- list()
edge_rows <- list()

for (slug in names(timepoint_lookup)) {
  timepoint <- unname(timepoint_lookup[[slug]])
  summary_path <- file.path(timepoint_root, slug, "fold_summary.tsv")
  summary <- fread(summary_path)
  if (nrow(summary) != expected_folds[[slug]] ||
    anyDuplicated(summary$heldout_replicate)) {
    stop("Unexpected fold contract at ", slug)
  }
  fold_rows[[slug]] <- summary

  per_fold <- vector("list", nrow(summary))
  for (fold_index in seq_len(nrow(summary))) {
    replicate <- summary$heldout_replicate[[fold_index]]
    decomposition_path <- file.path(
      timepoint_root, slug, "folds", replicate,
      "full_edge_decomposition.tsv.gz"
    )
    x <- fread(decomposition_path)
    if (!all(c(
      "regulator", "target", "full_weight", "full_rank",
      required_state_columns
    ) %in% names(x)) || anyDuplicated(x[, .(regulator, target)]) ||
      nrow(x) != summary$full_edges[[fold_index]]) {
      stop("Invalid decomposition file: ", decomposition_path)
    }
    for (column in required_state_columns) {
      if (anyNA(x[[column]]) || !is.logical(x[[column]])) {
        stop("Invalid logical state ", column, " in ", decomposition_path)
      }
    }
    if (any(x$fit_eligible &
      !(x$rna_model_eligible & x$multiomic_candidate_present)) ||
      any(x$selected_in_fold & !x$fit_eligible)) {
      stop("Non-nested decomposition states in ", decomposition_path)
    }
    x[, `:=`(
      timepoint = timepoint,
      heldout_replicate = replicate,
      candidate_and_rna_model =
        multiomic_candidate_present & rna_model_eligible
    )]
    per_fold[[fold_index]] <- x
  }
  timepoint_edges <- rbindlist(per_fold)
  edge_rows[[slug]] <- timepoint_edges[, .(
    full_weight = unique(full_weight),
    full_rank = unique(full_rank),
    folds = .N,
    rna_detection_eligible_folds = sum(rna_detection_eligible),
    rna_model_eligible_folds = sum(rna_model_eligible),
    candidate_present_folds = sum(multiomic_candidate_present),
    candidate_and_rna_model_folds = sum(candidate_and_rna_model),
    fit_eligible_folds = sum(fit_eligible),
    selected_folds = sum(selected_in_fold)
  ), by = .(timepoint, regulator, target)]
}

folds <- rbindlist(fold_rows)
edges <- rbindlist(edge_rows)
if (nrow(folds) != 14L || any(edges$selected_folds > edges$fit_eligible_folds) ||
  any(edges$fit_eligible_folds > edges$candidate_and_rna_model_folds) ||
  any(edges$candidate_and_rna_model_folds > edges$rna_model_eligible_folds)) {
  stop("Aggregated decomposition contract failed.")
}
edges[, `:=`(
  rna_model_eligibility_fraction = rna_model_eligible_folds / folds,
  candidate_retention_given_rna_model = safe_fraction(
    candidate_and_rna_model_folds, rna_model_eligible_folds
  ),
  fit_eligibility_given_rna_and_candidate = safe_fraction(
    fit_eligible_folds, candidate_and_rna_model_folds
  ),
  selection_given_fit_eligible = safe_fraction(
    selected_folds, fit_eligible_folds
  ),
  unconditional_selection_fraction = selected_folds / folds
)]

summarize_opportunities <- function(x, label, value) {
  data.table(
    stratum = label,
    stratum_value = as.character(value),
    edges = nrow(x),
    fold_opportunities = sum(x$folds),
    rna_detection_eligible = sum(x$rna_detection_eligible_folds),
    rna_model_eligible = sum(x$rna_model_eligible_folds),
    candidate_and_rna_model = sum(x$candidate_and_rna_model_folds),
    fit_eligible = sum(x$fit_eligible_folds),
    selected = sum(x$selected_folds),
    rna_detection_fraction = safe_fraction(
      sum(x$rna_detection_eligible_folds), sum(x$folds)
    ),
    rna_model_fraction = safe_fraction(
      sum(x$rna_model_eligible_folds), sum(x$folds)
    ),
    candidate_retention_given_rna_model = safe_fraction(
      sum(x$candidate_and_rna_model_folds), sum(x$rna_model_eligible_folds)
    ),
    fit_eligibility_given_rna_and_candidate = safe_fraction(
      sum(x$fit_eligible_folds), sum(x$candidate_and_rna_model_folds)
    ),
    selection_given_fit_eligible = safe_fraction(
      sum(x$selected_folds), sum(x$fit_eligible_folds)
    ),
    unconditional_selection_fraction = safe_fraction(
      sum(x$selected_folds), sum(x$folds)
    )
  )
}

top_k_requested <- c(100L, 1000L, 10000L, 50000L, 100000L)
top_rank <- edges[,
  {
    maximum_rank <- max(full_rank)
    rbindlist(lapply(unique(pmin(top_k_requested, maximum_rank)), function(k) {
      summarize_opportunities(.SD[full_rank <= k], "top_k", k)
    }))
  },
  by = timepoint
]

edges[, rank_decile := pmin(10L, ceiling(10 * full_rank / .N)), by = timepoint]
rank_decile <- edges[, summarize_opportunities(
  .SD, "rank_decile", unique(rank_decile)
), by = .(timepoint, rank_decile)]
rank_decile[, rank_decile := NULL]

layer_attrition <- rbindlist(lapply(seq_len(nrow(folds)), function(i) {
  x <- folds[i]
  counts <- c(
    complete_timepoint_edges = x$full_edges,
    rna_detection_eligible = x$full_edge_rna_detection_eligible,
    rna_model_eligible = x$full_edge_rna_model_eligible,
    multiomic_candidate_and_rna_model =
      x$full_edge_candidate_and_rna_model,
    fit_eligible = x$full_edge_fit_eligible,
    selected = x$full_edge_selected
  )
  data.table(
    dataset = x$dataset,
    timepoint = x$timepoint,
    heldout_replicate = x$heldout_replicate,
    layer = names(counts),
    edges = as.integer(counts),
    fraction_of_complete_edges = as.numeric(counts) / x$full_edges
  )
}))

fwrite(folds, file.path(stage_root, "all_fold_metrics.tsv"), sep = "\t")
fwrite(edges, file.path(
  stage_root, "full_edge_multilayer_stability.tsv.gz"
), sep = "\t")
fwrite(top_rank, file.path(stage_root, "top_rank_decomposition.tsv"), sep = "\t")
fwrite(rank_decile, file.path(
  stage_root, "rank_decile_decomposition.tsv"
), sep = "\t")
fwrite(layer_attrition, file.path(stage_root, "layer_attrition.tsv"), sep = "\t")

if (!file.rename(stage_root, output_root)) {
  stop("Unable to promote rebuilt-domain decomposition.")
}
cat("gse274113_rebuilt_domain_decomposition\tPASS\t", output_root, "\n", sep = "")
