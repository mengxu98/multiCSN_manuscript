#!/usr/bin/env Rscript

source("code/functions/table_io.R")

source("code/functions/runtime.R")
configure_task_runtime()


# ---- Parameters ----
# seeds: 1001:1010.
arguments <- list(
  input_root = "results/scmultisim/multiseed",
  output_root = "results/scmultisim/multiseed"
)
# ---------------------------------------------------------------

input_root <- normalizePath(arguments$input_root, mustWork = TRUE)
output_root <- normalizePath(arguments$output_root, mustWork = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

quantile_value <- function(x, p) as.numeric(stats::quantile(x, p, names = FALSE))

metrics <- utils::read.delim(file.path(input_root, "endpoint_metrics_all_seeds.tsv"))
seeds <- sort(unique(metrics$seed))
metrics$ranking_scheme <- "absolute_coefficient"
metrics$ap_over_prevalence <- metrics$average_precision / metrics$prevalence
score_files <- c(
  `TF-gene:static_input_GRN` = "scored_tf_gene_static.tsv.gz",
  `TF-gene:dynamic_union` = "scored_tf_gene_dynamic.tsv.gz",
  `region-gene:all_variable_regions` = "scored_region_gene.tsv.gz"
)
support <- do.call(rbind, lapply(seeds, function(seed) {
  do.call(rbind, lapply(names(score_files), function(scope) {
    scored <- utils::read.delim(file.path(
      input_root, paste0("seed_", seed), score_files[[scope]]
    ))
    selected <- scored$score > 0
    data.frame(
      seed = seed, scope = scope, selected = sum(selected),
      true_selected = sum(selected & scored$truth),
      selected_precision = sum(selected & scored$truth) / sum(selected),
      truth_recall = sum(selected & scored$truth) / sum(scored$truth),
      stringsAsFactors = FALSE
    )
  }))
}))
combined <- merge(metrics, support, by = c("seed", "scope"))
summary <- do.call(rbind, lapply(split(combined, combined$scope), function(x) {
  data.frame(
    scope = x$scope[1], seeds = nrow(x),
    prevalence_median = stats::median(x$prevalence),
    predicted_median = stats::median(x$n_predicted),
    true_selected_median = stats::median(x$true_selected),
    selected_precision_median = stats::median(x$selected_precision),
    truth_recall_median = stats::median(x$truth_recall),
    auroc_median = stats::median(x$auroc),
    auroc_q25 = quantile_value(x$auroc, 0.25),
    auroc_q75 = quantile_value(x$auroc, 0.75),
    average_precision_median = stats::median(x$average_precision),
    average_precision_q25 = quantile_value(x$average_precision, 0.25),
    average_precision_q75 = quantile_value(x$average_precision, 0.75),
    ap_over_prevalence_median = stats::median(x$ap_over_prevalence),
    stringsAsFactors = FALSE
  )
}))

static <- combined[combined$scope == "TF-gene:static_input_GRN", ]
dynamic <- combined[combined$scope == "TF-gene:dynamic_union", ]
static <- static[order(static$seed), ]
dynamic <- dynamic[order(dynamic$seed), ]
paired_truth <- data.frame(
  metric = c("average_precision", "auroc", "truth_recall"),
  median_static_minus_dynamic = c(
    stats::median(static$average_precision - dynamic$average_precision),
    stats::median(static$auroc - dynamic$auroc),
    stats::median(static$truth_recall - dynamic$truth_recall)
  ),
  stringsAsFactors = FALSE
)
write_tsv(combined, file.path(output_root, "per_seed_diagnostics.tsv"), writer = "table")
write_tsv(summary, file.path(output_root, "endpoint_diagnostic_summary.tsv"), writer = "table")
write_tsv(paired_truth, file.path(output_root, "static_dynamic_truth_sensitivity.tsv"), writer = "table")
cat("scMultiSim multiseed analysis: PASS\n")
