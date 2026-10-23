#!/usr/bin/env Rscript

source("code/functions/table_io.R")
summarize_pairwise_fold_stability <- function(timepoints = "results/gse274113/stability",
  output = "results/gse274113/stability") {

suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)

safe_spearman <- function(x, y) {
  if (length(x) < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(cor(x, y, method = "spearman"))
}


node_jaccard <- function(left, right, node) {
  left_degree <- left[, .(left_degree = .N), by = node]
  right_degree <- right[, .(right_degree = .N), by = node]
  shared <- merge(
    left[, .(regulator, target)], right[, .(regulator, target)],
    by = c("regulator", "target")
  )[, .(shared_degree = .N), by = node]
  result <- Reduce(function(a, b) merge(a, b, by = node, all = TRUE), list(
    left_degree, right_degree, shared
  ))
  for (column in c("left_degree", "right_degree", "shared_degree")) {
    set(result, which(is.na(result[[column]])), column, 0L)
  }
  result[, union_degree := left_degree + right_degree - shared_degree]
  result[union_degree > 0L, shared_degree / union_degree]
}

# ---- Parameters ----
# timepoints and output directory configuration.
arguments <- list(timepoints = timepoints, output = output)
# ---------------------------------------------------------------
timepoint_root <- normalizePath(arguments$timepoints, mustWork = TRUE)
final_root <- normalizePath(arguments$output, mustWork = FALSE)
dir.create(dirname(final_root), recursive = TRUE, showWarnings = FALSE)
output_root <- tempfile(".result-files-", tmpdir = dirname(final_root))
dir.create(output_root, recursive = FALSE)

timepoint_lookup <- c(
  day_7 = "day 7", day_9 = "day 9", day_11 = "day 11", day_14 = "day 14"
)
expected_folds <- c(day_7 = 4L, day_9 = 4L, day_11 = 3L, day_14 = 3L)
pair_rows <- list()
consensus_rows <- list()
row_index <- 0L
for (slug in names(timepoint_lookup)) {
  folds_root <- file.path(timepoint_root, "folds")
  fold_dirs <- sort(list.dirs(folds_root, full.names = TRUE, recursive = FALSE))
  fold_dirs <- fold_dirs[startsWith(basename(fold_dirs), paste0(slug, "_"))]
  if (length(fold_dirs) != expected_folds[[slug]]) {
    stop("Unexpected fold directory count at ", slug)
  }
  fold_names <- substring(basename(fold_dirs), nchar(slug) + 2L)
  networks <- setNames(lapply(fold_dirs, function(root) {
    x <- fread(file.path(root, "tf_gene.tsv"))[, .(
      regulator = as.character(regulator), target = as.character(target),
      weight = as.numeric(weight)
    )]
    if (!nrow(x) || anyDuplicated(x[, .(regulator, target)]) ||
      any(!is.finite(x$weight)) || any(x$weight == 0)) {
      stop("Invalid fold network: ", root)
    }
    x
  }), fold_names)

  pairs <- combn(fold_names, 2L, simplify = FALSE)
  for (pair in pairs) {
    left <- networks[[pair[[1]]]]
    right <- networks[[pair[[2]]]]
    shared <- merge(left, right,
      by = c("regulator", "target"),
      suffixes = c("_left", "_right")
    )
    target_jaccard <- node_jaccard(left, right, "target")
    regulator_jaccard <- node_jaccard(left, right, "regulator")
    row_index <- row_index + 1L
    pair_rows[[row_index]] <- data.frame(
      dataset = "GSE274113", timepoint = timepoint_lookup[[slug]],
      fold_left = pair[[1]], fold_right = pair[[2]],
      left_edges = nrow(left), right_edges = nrow(right),
      shared_edges = nrow(shared),
      edge_jaccard = nrow(shared) / (nrow(left) + nrow(right) - nrow(shared)),
      shared_edge_sign_agreement = mean(
        sign(shared$weight_left) == sign(shared$weight_right)
      ),
      shared_edge_weight_spearman = safe_spearman(
        abs(shared$weight_left), abs(shared$weight_right)
      ),
      median_target_jaccard = median(target_jaccard),
      mean_target_jaccard = mean(target_jaccard),
      median_regulator_jaccard = median(regulator_jaccard),
      mean_regulator_jaccard = mean(regulator_jaccard),
      stringsAsFactors = FALSE
    )
  }

  stacked <- rbindlist(lapply(names(networks), function(fold) {
    copy(networks[[fold]])[, fold := fold]
  }))
  consensus <- stacked[, .(
    selected_folds = uniqueN(fold),
    positive_folds = sum(weight > 0), negative_folds = sum(weight < 0),
    median_absolute_weight = median(abs(weight)),
    minimum_absolute_weight = min(abs(weight)),
    maximum_absolute_weight = max(abs(weight))
  ), by = .(regulator, target)]
  consensus[, `:=`(
    total_folds = length(networks),
    selection_fraction = selected_folds / length(networks),
    sign_consistent_when_selected =
      positive_folds == selected_folds | negative_folds == selected_folds
  )]
  setorder(consensus, -selected_folds, regulator, target)
  fwrite(consensus, file.path(
    output_root, paste0(slug, "_edge_fold_support.tsv.gz")
  ), sep = "\t")
  consensus_rows[[slug]] <- data.frame(
    dataset = "GSE274113", timepoint = timepoint_lookup[[slug]],
    folds = length(networks), union_edges = nrow(consensus),
    all_fold_edges = sum(consensus$selected_folds == length(networks)),
    all_fold_fraction_of_union = mean(
      consensus$selected_folds == length(networks)
    ),
    at_least_folds_minus_one_edges = sum(
      consensus$selected_folds >= length(networks) - 1L
    ),
    at_least_folds_minus_one_fraction_of_union = mean(
      consensus$selected_folds >= length(networks) - 1L
    ),
    selected_more_than_once_edges = sum(consensus$selected_folds >= 2L),
    selected_more_than_once_fraction_of_union = mean(
      consensus$selected_folds >= 2L
    ),
    sign_consistent_all_fold_edges = mean(
      consensus$sign_consistent_when_selected[
        consensus$selected_folds == length(networks)
      ]
    ), stringsAsFactors = FALSE
  )
}

pairs <- rbindlist(pair_rows)
consensus_summary <- rbindlist(consensus_rows)
metrics <- c(
  "edge_jaccard", "shared_edge_sign_agreement", "shared_edge_weight_spearman",
  "median_target_jaccard", "mean_target_jaccard",
  "median_regulator_jaccard", "mean_regulator_jaccard"
)
if (nrow(pairs) != 18L || any(!vapply(pairs[, ..metrics], function(x) {
  all(is.finite(x)) && all(x >= 0 & x <= 1)
}, logical(1)))) {
  stop("Pairwise fold stability contract failed.")
}
pair_summary <- melt(
  pairs,
  id.vars = "timepoint", measure.vars = metrics,
  variable.name = "metric", value.name = "value"
)[, .(
  fold_pairs = .N, minimum = min(value), first_quartile = quantile(value, 0.25),
  median = median(value), third_quartile = quantile(value, 0.75),
  maximum = max(value), mean = mean(value), standard_deviation = sd(value)
), by = .(timepoint, metric)]

fwrite(pairs, file.path(output_root, "pairwise_fold_metrics.tsv"), sep = "\t")
fwrite(pair_summary, file.path(output_root, "pairwise_fold_summary.tsv"), sep = "\t")
fwrite(consensus_summary, file.path(output_root, "edge_fold_support_summary.tsv"),
  sep = "\t"
)
publish_result_files(output_root, final_root)
cat("gse274113_pairwise_fold_stability\tPASS\t", final_root, "\n", sep = "")
}



if (sys.nframe() == 0L) summarize_pairwise_fold_stability()
