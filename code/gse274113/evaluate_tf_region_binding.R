#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})
data.table::setDTthreads(1L)


regions_to_granges <- function(values) {
  values <- as.character(values)
  matched <- regexec("^(.+?)[-:]([0-9]+)-([0-9]+)$", values, perl = TRUE)
  pieces <- regmatches(values, matched)
  valid <- lengths(pieces) == 4L
  if (!all(valid)) stop("Unparseable region coordinates.")
  chrom <- vapply(pieces, `[[`, character(1), 2L)
  start <- as.integer(vapply(pieces, `[[`, character(1), 3L))
  end <- as.integer(vapply(pieces, `[[`, character(1), 4L))
  GenomicRanges::GRanges(chrom, IRanges::IRanges(start = start, end = end))
}

table_to_granges <- function(table) {
  GenomicRanges::GRanges(
    as.character(table$chrom),
    IRanges::IRanges(start = as.integer(table$start), end = as.integer(table$end))
  )
}

set_metrics <- function(labels, predicted) {
  labels <- as.logical(labels)
  predicted <- as.logical(predicted)
  tp <- sum(labels & predicted)
  fp <- sum(!labels & predicted)
  fn <- sum(labels & !predicted)
  tn <- sum(!labels & !predicted)
  precision <- if (tp + fp) tp / (tp + fp) else NA_real_
  recall <- if (tp + fn) tp / (tp + fn) else NA_real_
  f1 <- if (is.finite(precision) && is.finite(recall) && precision + recall > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }
  prevalence <- mean(labels)
  fisher_p <- if ((tp + fp) && (tp + fn) && (tn + fp) && (tn + fn)) {
    stats::fisher.test(
      matrix(c(tp, fp, fn, tn), nrow = 2),
      alternative = "greater"
    )$p.value
  } else {
    NA_real_
  }
  data.frame(
    evaluable_regions = length(labels),
    evidence_positive_regions = tp + fn,
    predicted_regions_in_universe = tp + fp,
    true_positive = tp,
    false_positive = fp,
    false_negative = fn,
    true_negative = tn,
    precision = precision,
    recall = recall,
    f1 = f1,
    jaccard = if (tp + fp + fn) tp / (tp + fp + fn) else NA_real_,
    enrichment_over_prevalence = if (is.finite(precision) && prevalence > 0) {
      precision / prevalence
    } else {
      NA_real_
    },
    fisher_exact_greater_p = fisher_p,
    stringsAsFactors = FALSE
  )
}

# ---- Parameters ----
# External evidence was rebuilt on 2026-09-14; source URLs, hashes and conversion
# records: data/gse274113/raw_external_downloads/README.md.
# The four-method reproduction runner uses recovered server-131 LINGER inputs.
evaluate_tf_region_binding <- function(arguments = list(
  endpoints = "results/gse274113/networks",
  evidence = "data/gse274113/external_resources",
  truth = "data/gse274113/heldout_pseudobulk_full",
  output = "results/gse274113/multiomic_endpoint_evidence/TF_region_external_binding",
  methods = NULL
)) {
endpoint_root <- normalizePath(arguments$endpoints, mustWork = TRUE)
evidence_root <- normalizePath(arguments$evidence, mustWork = TRUE)
truth_root <- normalizePath(arguments$truth, mustWork = TRUE)
output_root <- normalizePath(arguments$output, mustWork = FALSE)
if (file.exists(output_root)) {
  stop("Output already exists and will not be overwritten: ", output_root)
}
dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
stage_root <- paste0(
  output_root, ".staging_", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
)
dir.create(stage_root, recursive = FALSE)

chip_path <- file.path(evidence_root, "primary_erythroid_ChIP_peaks_hg38.tsv.gz")
unibind_path <- file.path(evidence_root, "UniBind_perturbed_TFs_hg38.tsv.gz")
effect_summary_path <- file.path(truth_root, "effect_model_summary.tsv")
if (!all(file.exists(c(chip_path, unibind_path, effect_summary_path)))) {
  stop("Prepared binding evidence and held-out truth are required.")
}
chip <- fread(chip_path)
unibind <- fread(unibind_path)
required_binding <- c("regulator", "chrom", "start", "end")
if (!all(required_binding %in% names(chip)) ||
  !all(required_binding %in% names(unibind))) {
  stop("Prepared binding-evidence schema changed.")
}
effect_summary <- fread(effect_summary_path)
targets <- sort(unique(as.character(effect_summary$target)))
timepoints <- c("day 7", "day 9", "day 11", "day 14")
methods <- if (is.null(arguments$methods)) {
  c("multiCSN", "SCENICplus", "Pando")
} else {
  strsplit(arguments$methods, ",", fixed = TRUE)[[1]]
}
if (!length(methods) || anyDuplicated(methods) || any(!nzchar(methods))) {
  stop("methods must contain unique method names.")
}
if (length(targets) != 19L) stop("Expected 19 perturbed TFs.")

evidence_sets <- list(
  primary_erythroid_ChIP = chip,
  UniBind_direct_TFBS = unibind
)
evidence_ranges <- lapply(evidence_sets, function(table) {
  split(table, table$regulator) |>
    lapply(function(subset) GenomicRanges::reduce(table_to_granges(subset)))
})

rows <- list()
row_index <- 0L
for (timepoint in timepoints) {
  slug <- gsub(" ", "_", timepoint, fixed = TRUE)
  predictions <- list()
  for (method in methods) {
    path <- file.path(endpoint_root, method, slug, "tf_region.tsv.gz")
    if (!file.exists(path)) stop("Missing TF-region endpoint: ", path)
    predictions[[method]] <- fread(path)
  }
  for (regulator in targets) {
    regulator_id <- regulator
    truth_path <- file.path(
      truth_root, "atac_effects", paste0(slug, "__", regulator, ".tsv.gz")
    )
    if (!file.exists(truth_path)) stop("Missing ATAC universe: ", truth_path)
    universe <- fread(truth_path, select = "feature")$feature
    if (anyDuplicated(universe)) stop("Duplicated ATAC universe regions.")
    universe_ranges <- regions_to_granges(universe)
    for (evidence_name in names(evidence_ranges)) {
      if (!regulator %in% names(evidence_ranges[[evidence_name]])) next
      labels <- GenomicRanges::countOverlaps(
        universe_ranges, evidence_ranges[[evidence_name]][[regulator]],
        ignore.strand = TRUE
      ) > 0
      for (method in methods) {
        prediction <- predictions[[method]][regulator == regulator_id]
        predicted_total <- nrow(prediction)
        if (predicted_total) {
          predicted_ranges <- regions_to_granges(prediction$region)
          predicted <- GenomicRanges::countOverlaps(
            universe_ranges, predicted_ranges,
            ignore.strand = TRUE
          ) > 0
        } else {
          predicted <- rep(FALSE, length(universe_ranges))
        }
        row_index <- row_index + 1L
        metric <- set_metrics(labels, predicted)
        metric$method <- method
        metric$timepoint <- timepoint
        metric$regulator <- regulator
        metric$evidence_source <- evidence_name
        metric$predicted_relations_total <- predicted_total
        metric$score_policy <- "binary_set_only"
        rows[[row_index]] <- metric
      }
    }
  }
}

metrics <- rbindlist(rows, use.names = TRUE, fill = TRUE)
if (!nrow(metrics) || any(metrics$score_policy != "binary_set_only")) {
  stop("External TF-region evidence evaluation is empty or mis-scored.")
}
metrics[, `:=`(
  fisher_bh_q_endpoint_family = p.adjust(
    fisher_exact_greater_p,
    method = "BH"
  ),
  fisher_holm_p_endpoint_family = p.adjust(
    fisher_exact_greater_p,
    method = "holm"
  )
), by = evidence_source]
setcolorder(metrics, c(
  "method", "timepoint", "regulator", "evidence_source", "score_policy",
  setdiff(names(metrics), c(
    "method", "timepoint", "regulator", "evidence_source", "score_policy"
  ))
))
write_tsv(metrics, file.path(stage_root, "TF_region_external_binding_metrics.tsv.gz"))
metric_columns <- c(
  "precision", "recall", "f1", "jaccard", "enrichment_over_prevalence"
)
summarize_binary <- function(table, common = FALSE) {
  table[,
    {
      means <- setNames(lapply(.SD, function(x) {
        x <- x[is.finite(x)]
        if (length(x)) mean(x) else NA_real_
      }), paste0("macro_", metric_columns))
      finite_counts <- setNames(
        lapply(.SD, function(x) sum(is.finite(x))),
        paste0(metric_columns, "_finite_strata")
      )
      total_tp <- sum(true_positive)
      total_predicted <- sum(predicted_regions_in_universe)
      total_evidence <- sum(evidence_positive_regions)
      micro_precision <- if (total_predicted) total_tp / total_predicted else NA_real_
      micro_recall <- if (total_evidence) total_tp / total_evidence else NA_real_
      c(list(
        strata_total = .N,
        nonempty_prediction_strata = sum(predicted_regions_in_universe > 0L),
        true_positive_strata = sum(true_positive > 0L),
        predicted_regions_in_universe = total_predicted,
        evidence_positive_regions = total_evidence,
        true_positive = total_tp,
        micro_precision = micro_precision,
        micro_recall = micro_recall,
        micro_f1 = if (is.finite(micro_precision) && is.finite(micro_recall) &&
          micro_precision + micro_recall > 0) {
          2 * micro_precision * micro_recall / (micro_precision + micro_recall)
        } else {
          NA_real_
        },
        BH_significant_strata = sum(
          fisher_bh_q_endpoint_family < 0.05,
          na.rm = TRUE
        ),
        Holm_significant_strata = sum(
          fisher_holm_p_endpoint_family < 0.05,
          na.rm = TRUE
        )
      ), means, finite_counts)
    },
    by = .(method, evidence_source),
    .SDcols = metric_columns
  ]
}
macro <- summarize_binary(metrics)
write_tsv(macro, file.path(stage_root, "TF_region_external_binding_macro.tsv"))

common_keys <- metrics[, .(
  methods_observed = uniqueN(method),
  methods_nonempty = sum(predicted_regions_in_universe > 0L)
), by = .(timepoint, regulator, evidence_source)][
  methods_observed == length(methods) & methods_nonempty == length(methods),
  .(timepoint, regulator, evidence_source)
]
common_metrics <- merge(
  metrics, common_keys,
  by = c("timepoint", "regulator", "evidence_source"),
  all = FALSE
)
common_macro <- summarize_binary(common_metrics)
setnames(common_macro, "strata_total", "common_nonempty_strata")
write_tsv(
  common_macro,
  file.path(stage_root, "TF_region_external_binding_common_nonempty_macro.tsv")
)

if (!file.rename(stage_root, output_root)) {
  stop("Unable to atomically publish TF-region external evaluation.")
}
cat("TF_region_external_binding_evaluation\tPASS\t", output_root, "\n", sep = "")

}

if (sys.nframe() == 0L) evaluate_tf_region_binding()
