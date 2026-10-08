#!/usr/bin/env Rscript

source("code/functions/ranking_metrics.R")

suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)

# ---- Parameters ----
# Parameters: NETWORK_ROOT, TRUTH_ROOT, FORMAL_METRICS, OUTPUT_ROOT.
arguments <- list(
  network_root = "runs/gse274113/LINGER",
  truth_root = "data/gse274113/heldout_pseudobulk_full",
  formal_metrics = "results/gse274113/method_comparison/stratum_metrics.tsv",
  output = "results/gse274113/method_comparison"
)
# ---------------------------------------------------------------
network_root <- normalizePath(arguments$network_root, mustWork = TRUE)
truth_root <- normalizePath(arguments$truth_root, mustWork = TRUE)
formal_metrics_path <- normalizePath(arguments$formal_metrics, mustWork = TRUE)
result_root <- normalizePath(arguments$output, mustWork = FALSE)
timepoints <- c("day 7", "day 9", "day 11", "day 14")

formal <- fread(formal_metrics_path)
targets <- sort(unique(formal[endpoint == "TF-gene"]$regulator))
if (length(targets) != 19L) stop("Expected 19 perturbed TFs in the formal comparison.")

rows <- list()
cursor <- 1L
for (timepoint in timepoints) {
  slug <- gsub(" ", "_", timepoint, fixed = TRUE)
  network_path <- file.path(network_root, slug, "networks", "tf_gene.tsv.gz")
  if (!file.exists(network_path)) stop("Missing LINGER network: ", network_path)
  network <- fread(network_path)
  if (!all(c("regulator", "target", "score") %in% names(network))) {
    stop("LINGER TF-gene schema changed at ", slug)
  }
  network <- network[, .(
    regulator = as.character(regulator),
    feature = as.character(target),
    score = as.numeric(score)
  )]
  if (any(!is.finite(network$score)) || any(network$score <= 0)) {
    stop("LINGER exported a non-finite or non-positive TF-gene score at ", slug)
  }
  exported_rows <- nrow(network)
  network <- network[, .(score = max(score)), by = .(regulator, feature)]
  duplicate_pair_rows <- exported_rows - nrow(network)
  if (duplicate_pair_rows) stop("LINGER exported duplicate TF-gene pairs at ", slug)

  for (tf in targets) {
    truth_path <- file.path(truth_root, "rna_effects", paste0(slug, "__", tf, ".tsv.gz"))
    truth <- fread(truth_path)
    if (!all(c("feature", "logFC", "FDR") %in% names(truth)) || anyDuplicated(truth$feature)) {
      stop("Held-out RNA truth schema changed for ", tf, " at ", timepoint)
    }
    prediction <- network[regulator == tf, .(feature, score)]
    joined <- merge(
      truth[, .(feature = as.character(feature), FDR = as.numeric(FDR))],
      prediction,
      by = "feature", all.x = TRUE, sort = FALSE
    )[feature != tf & is.finite(FDR)]
    joined[, relation_present := !is.na(score)]
    joined[is.na(score), score := 0]
    joined[, truth_response := FDR <= 0.05]
    positives <- sum(joined$truth_response)
    prevalence <- positives / nrow(joined)
    ap <- threshold_average_precision(joined$truth_response, joined$score)
    selected_truth <- joined$truth_response & joined$score > 0
    rows[[cursor]] <- data.table(
      method = "LINGER",
      endpoint = "TF-gene",
      comparison_universe = "full_truth_feature_universe",
      regulator = tf,
      timepoint = timepoint,
      evaluable_features = nrow(joined),
      truth_responses_fdr_0_05 = positives,
      truth_response_prevalence = prevalence,
      native_relations = sum(joined$relation_present),
      score_available_relations = sum(joined$relation_present),
      predicted_nonzero = sum(joined$score > 0),
      direction_known_nonzero = 0L,
      prediction_density = mean(joined$score > 0),
      truth_response_coverage = if (positives) sum(selected_truth) / positives else NA_real_,
      response_average_precision = ap,
      response_ap_lift_over_prevalence = if (prevalence > 0) ap / prevalence else NA_real_,
      response_auroc = tie_safe_auroc(joined$truth_response, joined$score),
      response_precision_at_truth_k = expected_precision_at_k(
        joined$truth_response, joined$score, positives
      ),
      directional_precision_at_truth_k = NA_real_,
      direction_known_coverage_selected_truth = NA_real_,
      direction_accuracy_selected_truth = NA_real_,
      signed_spearman_all = NA_real_,
      signed_spearman_selected = NA_real_
    )
    cursor <- cursor + 1L
  }
}

linger <- rbindlist(rows)
if (nrow(linger) != 76L || anyDuplicated(linger[, .(regulator, timepoint)]) ||
    any(linger$predicted_nonzero > linger$evaluable_features)) {
  stop("LINGER evaluation did not produce the declared 19 TF x 4 time-point design.")
}
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
fwrite(linger, file.path(result_root, "linger_stratum_metrics.tsv"), sep = "\t", na = "NA")
combined <- rbindlist(list(formal[endpoint == "TF-gene"], linger), use.names = TRUE)
fwrite(combined, file.path(result_root, "combined_tf_gene_stratum_metrics.tsv"), sep = "\t", na = "NA")

metric_labels <- c(
  response_average_precision = "Response AP",
  response_auroc = "AUROC",
  truth_response_coverage = "Truth-response coverage",
  prediction_density = "Prediction density"
)
base <- formal[method == "multiCSN" & endpoint == "TF-gene"]
paired <- merge(
  base,
  linger,
  by = c("endpoint", "comparison_universe", "regulator", "timepoint"),
  suffixes = c("_multicsn", "_linger"), sort = TRUE
)
if (nrow(paired) != 76L) stop("The LINGER and multiCSN strata are not paired.")
source_rows <- rbindlist(lapply(names(metric_labels), function(metric) {
  data.table(
    regulator = paired$regulator,
    timepoint = paired$timepoint,
    metric = metric,
    metric_label = unname(metric_labels[[metric]]),
    multicsn = paired[[paste0(metric, "_multicsn")]],
    linger = paired[[paste0(metric, "_linger")]],
    unit = "TF x time point"
  )
}))
tf_means <- source_rows[
  is.finite(multicsn) & is.finite(linger),
  .(multicsn = mean(multicsn), linger = mean(linger)),
  by = .(regulator, metric, metric_label)
][, `:=`(
  timepoint = "Mean of finite time points",
  unit = "TF mean"
)]
figure_source <- rbindlist(list(source_rows, tf_means), use.names = TRUE)
figure_source[, difference_multicsn_minus_linger := multicsn - linger]
fwrite(figure_source, file.path(result_root, "figure_source_data.tsv"), sep = "\t", na = "NA")
cat("linger_tf_gene_evaluation\tPASS\t", result_root, "\n", sep = "")
