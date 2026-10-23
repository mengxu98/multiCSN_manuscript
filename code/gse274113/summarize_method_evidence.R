#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)


# ---- Parameters ----
# input: stratum_metrics.tsv directory; output: method evidence summary directory.
arguments <- list(
  input = "results/gse274113/comparison",
  output = "results/gse274113/comparison"
)
# ---------------------------------------------------------------
input_root <- normalizePath(arguments$input, mustWork = TRUE)
output_root <- normalizePath(arguments$output, mustWork = FALSE)
dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
stage_root <- tempfile(".result-files-", tmpdir = dirname(output_root))
dir.create(stage_root, recursive = FALSE)

metric_path <- file.path(input_root, "stratum_metrics.tsv")
metrics <- fread(metric_path)

required <- c(
  "method", "endpoint", "comparison_universe", "regulator", "timepoint",
  "truth_response_coverage", "response_average_precision", "response_auroc",
  "response_ap_lift_over_prevalence", "response_precision_at_truth_k",
  "directional_precision_at_truth_k", "direction_accuracy_selected_truth",
  "signed_spearman_all"
)
if (!all(required %in% names(metrics)) ||
  anyDuplicated(metrics[, .(
    method, endpoint, comparison_universe,
    regulator, timepoint
  )])) {
  stop("The registered stratum table is incomplete or contains duplicate keys.")
}

metric_columns <- c(
  "truth_response_coverage", "response_average_precision", "response_auroc",
  "response_ap_lift_over_prevalence", "response_precision_at_truth_k",
  "directional_precision_at_truth_k", "direction_accuracy_selected_truth",
  "signed_spearman_all"
)
comparators <- c("Pando", "SCENIC+")
universes <- "full_truth_feature_universe"

difference_rows <- list()
test_rows <- list()
cursor <- 0L
test_cursor <- 0L
for (universe in universes) {
  reference <- metrics[
    method == "multiCSN" & endpoint == "TF-gene" &
      comparison_universe == universe
  ]
  if (nrow(reference) != 76L) stop("Incomplete multiCSN strata in ", universe)
  for (comparator in comparators) {
    other <- metrics[
      method == comparator & endpoint == "TF-gene" &
        comparison_universe == universe
    ]
    paired <- merge(
      reference, other,
      by = c("endpoint", "comparison_universe", "regulator", "timepoint"),
      suffixes = c("_multicsn", "_comparator"), sort = TRUE
    )
    if (nrow(paired) != 76L) stop("Incomplete paired strata for ", comparator)
    for (metric in metric_columns) {
      left <- paired[[paste0(metric, "_multicsn")]]
      right <- paired[[paste0(metric, "_comparator")]]
      paired[, difference__ := right - left]
      cluster <- paired[is.finite(difference__), .(
        finite_timepoints = .N,
        mean_difference = mean(difference__),
        median_difference = median(difference__)
      ), by = regulator]
      cluster <- merge(
        data.table(regulator = sort(unique(reference$regulator))),
        cluster,
        by = "regulator", all.x = TRUE, sort = TRUE
      )
      cursor <- cursor + 1L
      difference_rows[[cursor]] <- cluster[, .(
        comparator,
        reference = "multiCSN", endpoint = "TF-gene",
        comparison_universe = universe, metric,
        difference_definition = "comparator minus multiCSN",
        regulator, finite_timepoints, mean_difference, median_difference
      )]
      finite_cluster <- cluster[is.finite(mean_difference)]
      difference <- finite_cluster$mean_difference
      wilcoxon_p <- if (length(difference) >= 3L && any(difference != 0)) {
        suppressWarnings(wilcox.test(
          difference,
          mu = 0, paired = FALSE, exact = FALSE,
          conf.int = FALSE
        )$p.value)
      } else {
        NA_real_
      }
      test_cursor <- test_cursor + 1L
      test_rows[[test_cursor]] <- data.table(
        comparator,
        reference = "multiCSN", endpoint = "TF-gene",
        comparison_universe = universe, metric,
        difference_definition = "comparator minus multiCSN",
        independent_unit = "perturbed TF mean across its finite timepoints",
        regulator_clusters_total = 19L,
        regulator_clusters_finite = nrow(finite_cluster),
        mean_cluster_difference = if (length(difference)) mean(difference) else NA_real_,
        median_cluster_difference = if (length(difference)) median(difference) else NA_real_,
        comparator_wins = sum(difference > 0),
        exact_ties = sum(difference == 0),
        comparator_losses = sum(difference < 0),
        wilcoxon_signed_rank_p = wilcoxon_p
      )
    }
  }
}

cluster_differences <- rbindlist(difference_rows)
tests <- rbindlist(test_rows)
tests[, holm_p_two_comparators := p.adjust(
  wilcoxon_signed_rank_p,
  method = "holm"
), by = .(endpoint, comparison_universe, metric)]
tests[, multiplicity_family :=
  "two comparator-versus-multiCSN tests within endpoint, universe and metric"]

if (nrow(cluster_differences) != 2L * 8L * 19L ||
  nrow(tests) != 2L * 8L ||
  any(tests$regulator_clusters_finite > tests$regulator_clusters_total) ||
  any(tests$comparator_wins + tests$exact_ties + tests$comparator_losses !=
    tests$regulator_clusters_finite)) {
  stop("TF-cluster comparison accounting failed.")
}

write_tsv(cluster_differences, file.path(stage_root, "tf_cluster_metric_differences.tsv"))
write_tsv(tests, file.path(stage_root, "paired_method_cluster_tests.tsv"))

publish_result_files(stage_root, output_root)
cat("gse274113_current_method_evidence\tPASS\t", output_root, "\n", sep = "")
