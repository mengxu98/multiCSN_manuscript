#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)


# ---- Parameters ----
arguments <- list(
  full_timepoints = "results/gse274113/stability",
  fold_support = "results/gse274113/stability",
  output = "results/gse274113/stability"
)
# ---------------------------------------------------------------
full_root <- normalizePath(arguments$full_timepoints, mustWork = TRUE)
support_root <- normalizePath(arguments$fold_support, mustWork = TRUE)
final_root <- normalizePath(arguments$output, mustWork = FALSE)
dir.create(dirname(final_root), recursive = TRUE, showWarnings = FALSE)
stage_root <- tempfile(".result-files-", tmpdir = dirname(final_root))
dir.create(stage_root, recursive = FALSE)

timepoint_lookup <- c(
  day_7 = "day 7", day_9 = "day 9", day_11 = "day 11", day_14 = "day 14"
)
expected_folds <- c(day_7 = 4L, day_9 = 4L, day_11 = 3L, day_14 = 3L)
top_k_requested <- c(100L, 1000L, 10000L, 50000L, 100000L)
top_rows <- list()
decile_rows <- list()
edge_rows <- list()
for (slug in names(timepoint_lookup)) {
  full_path <- file.path(full_root, "full", slug, "tf_gene.tsv")
  support_path <- file.path(
    support_root, paste0(slug, "_edge_fold_support.tsv.gz")
  )
  full <- fread(full_path)[, .(
    regulator = as.character(regulator), target = as.character(target),
    weight = as.numeric(weight),
    deletion_delta_bic = as.numeric(deletion_delta_bic),
    target_detected_control_cells = as.integer(target_detected_control_cells),
    regulator_detected_control_cells = as.integer(regulator_detected_control_cells)
  )]
  support <- fread(support_path)[, .(
    regulator, target, selected_folds, total_folds, positive_folds,
    negative_folds, median_absolute_weight
  )]
  if (anyDuplicated(full[, .(regulator, target)]) ||
    anyDuplicated(support[, .(regulator, target)])) {
    stop("Duplicate edge in calibration input: ", slug)
  }
  edge <- merge(full, support, by = c("regulator", "target"), all.x = TRUE)
  if (nrow(edge) != nrow(full)) stop("Full network join changed row count.")
  edge[is.na(selected_folds), `:=`(
    selected_folds = 0L, total_folds = expected_folds[[slug]],
    positive_folds = 0L, negative_folds = 0L
  )]
  if (any(edge$total_folds != expected_folds[[slug]]) ||
    any(edge$selected_folds < 0L | edge$selected_folds > edge$total_folds)) {
    stop("Fold-support contract failed: ", slug)
  }
  edge[, weight_magnitude__ := abs(weight)]
  setorder(edge, -weight_magnitude__, regulator, target)
  edge[, weight_magnitude__ := NULL]
  edge[, full_rank := .I]
  edge[, full_rank_fraction := full_rank / .N]
  edge[, rank_decile := pmin(10L, ceiling(10 * full_rank_fraction))]
  edge[, selection_fraction := selected_folds / total_folds]
  edge[, sign_consistent_with_full := fifelse(
    selected_folds == 0L, NA,
    fifelse(weight > 0, negative_folds == 0L, positive_folds == 0L)
  )]
  edge[, timepoint := timepoint_lookup[[slug]]]

  available_top_k <- unique(pmin(top_k_requested, nrow(edge)))
  top_rows[[slug]] <- rbindlist(lapply(available_top_k, function(k) {
    x <- edge[seq_len(k)]
    data.frame(
      dataset = "GSE274113", timepoint = timepoint_lookup[[slug]],
      top_k = k, full_edges = nrow(edge), folds = expected_folds[[slug]],
      mean_selection_fraction = mean(x$selection_fraction),
      median_selection_fraction = median(x$selection_fraction),
      selected_in_any_fold_fraction = mean(x$selected_folds >= 1L),
      selected_in_at_least_two_folds_fraction = mean(x$selected_folds >= 2L),
      selected_in_all_folds_fraction = mean(
        x$selected_folds == expected_folds[[slug]]
      ),
      sign_consistency_with_full_among_reselected = mean(
        x$sign_consistent_with_full,
        na.rm = TRUE
      ), stringsAsFactors = FALSE
    )
  }))
  decile_rows[[slug]] <- edge[, .(
    edges = .N,
    minimum_rank = min(full_rank), maximum_rank = max(full_rank),
    mean_selection_fraction = mean(selection_fraction),
    median_selection_fraction = median(selection_fraction),
    selected_in_any_fold_fraction = mean(selected_folds >= 1L),
    selected_in_at_least_two_folds_fraction = mean(selected_folds >= 2L),
    selected_in_all_folds_fraction = mean(selected_folds == total_folds),
    sign_consistency_with_full_among_reselected = mean(
      sign_consistent_with_full,
      na.rm = TRUE
    ),
    median_target_detection = as.numeric(median(target_detected_control_cells)),
    median_regulator_detection = as.numeric(median(regulator_detected_control_cells))
  ), by = .(timepoint, rank_decile)]
  edge_rows[[slug]] <- edge[, .(
    timepoint, regulator, target, weight, deletion_delta_bic,
    target_detected_control_cells, regulator_detected_control_cells,
    full_rank, full_rank_fraction, rank_decile, selected_folds, total_folds,
    selection_fraction, sign_consistent_with_full, median_absolute_weight
  )]
}

top <- rbindlist(top_rows)
deciles <- rbindlist(decile_rows)
edges <- rbindlist(edge_rows)
if (any(!is.finite(top$mean_selection_fraction)) ||
  any(!is.finite(deciles$mean_selection_fraction)) ||
  any(edges$selection_fraction < 0 | edges$selection_fraction > 1)) {
  stop("Calibration output contract failed.")
}
fwrite(top, file.path(stage_root, "top_rank_stability.tsv"), sep = "\t")
fwrite(deciles, file.path(stage_root, "rank_decile_stability.tsv"), sep = "\t")
fwrite(edges, file.path(stage_root, "full_edge_fold_support.tsv.gz"), sep = "\t")
publish_result_files(stage_root, final_root)
cat("gse274113_weight_stability_calibration\tPASS\t", final_root, "\n", sep = "")
