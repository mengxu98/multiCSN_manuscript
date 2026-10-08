#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(multiCSN)
})
data.table::setDTthreads(1L)


# ---- Parameters ----
# input: complete Seurat object for NT controls.
# full_timepoints and folds inferred from replicate_deletion_stability directory structure.
arguments <- list(
  input = "data/gse274113/objects/GSE274113_control_full/GSE274113_all_NT_controls_full_seurat.rds",
  full_timepoints = "results/gse274113/replicate_deletion_stability/timepoints",
  folds = "results/gse274113/replicate_deletion_stability/timepoints",
  output = "results/gse274113/replicate_deletion_stability/eligibility_conditioned_stability"
)
# ---------------------------------------------------------------
input_path <- normalizePath(arguments$input, mustWork = TRUE)
full_root <- normalizePath(arguments$full_timepoints, mustWork = TRUE)
fold_root <- normalizePath(arguments$folds, mustWork = TRUE)
final_root <- normalizePath(arguments$output, mustWork = FALSE)
if (file.exists(final_root)) stop("Output already exists: ", final_root)
dir.create(dirname(final_root), recursive = TRUE, showWarnings = FALSE)
stage_root <- paste0(final_root, ".staging_", format(
  Sys.time(), "%Y%m%dT%H%M%SZ",
  tz = "UTC"
))
dir.create(stage_root, recursive = FALSE)

object <- readRDS(input_path)
params <- multiCSN::Params(object)
counts <- multiCSN::get_layer_data(
  object,
  assay = params$rna_assay, layer = "counts"
)
metadata <- object[[]]
timepoint_lookup <- c(
  day_7 = "day 7", day_9 = "day 9", day_11 = "day 11", day_14 = "day 14"
)
expected_folds <- c(day_7 = 4L, day_9 = 4L, day_11 = 3L, day_14 = 3L)
top_k_requested <- c(100L, 1000L, 10000L, 50000L, 100000L)
edge_rows <- list()
top_rows <- list()
decile_rows <- list()

for (slug in names(timepoint_lookup)) {
  timepoint <- timepoint_lookup[[slug]]
  full_path <- file.path(full_root, slug, "detected_20", "tf_gene.tsv")
  full <- fread(full_path)[, .(
    regulator = as.character(regulator), target = as.character(target),
    weight = as.numeric(weight),
    target_detected_control_cells = as.integer(target_detected_control_cells),
    regulator_detected_control_cells = as.integer(regulator_detected_control_cells)
  )]
  if (anyDuplicated(full[, .(regulator, target)])) stop("Duplicate full edge.")
  full[, weight_magnitude__ := abs(weight)]
  setorder(full, -weight_magnitude__, regulator, target)
  full[, `:=`(
    weight_magnitude__ = NULL, full_rank = .I,
    rank_decile = pmin(10L, ceiling(10 * .I / .N))
  )]
  cells <- rownames(metadata)[as.character(metadata$Timepoint) == timepoint]
  replicates <- sort(unique(as.character(metadata[cells, "replicate"])))
  if (length(replicates) != expected_folds[[slug]]) stop("Fold count mismatch.")
  genes <- sort(unique(c(full$regulator, full$target)))
  if (!all(genes %in% rownames(counts))) stop("Full edge gene absent from counts.")
  full_detection <- as.integer(Matrix::rowSums(counts[genes, cells] != 0))
  names(full_detection) <- genes
  if (!identical(
    unname(full_detection[full$target]),
    full$target_detected_control_cells
  ) ||
    !identical(
      unname(full_detection[full$regulator]),
      full$regulator_detected_control_cells
    )) {
    stop("Recorded full-data detection differs from source counts: ", slug)
  }

  selected_matrix <- matrix(FALSE, nrow(full), length(replicates),
    dimnames = list(NULL, replicates)
  )
  eligible_matrix <- matrix(FALSE, nrow(full), length(replicates),
    dimnames = list(NULL, replicates)
  )
  full_keys <- paste(full$regulator, full$target, sep = "\r")
  for (fold_index in seq_along(replicates)) {
    replicate <- replicates[[fold_index]]
    heldout_cells <- cells[as.character(metadata[cells, "replicate"]) == replicate]
    heldout_detection <- as.integer(Matrix::rowSums(
      counts[genes, heldout_cells] != 0
    ))
    names(heldout_detection) <- genes
    remaining_detection <- full_detection - heldout_detection
    eligible_matrix[, fold_index] <-
      remaining_detection[full$regulator] >= 20L &
        remaining_detection[full$target] >= 20L
    fold_path <- file.path(
      fold_root, slug, "folds", replicate, "tf_gene.tsv"
    )
    fold <- fread(fold_path, select = c("regulator", "target"))
    selected_matrix[, fold_index] <- full_keys %chin%
      paste(fold$regulator, fold$target, sep = "\r")
    if (any(selected_matrix[, fold_index] & !eligible_matrix[, fold_index])) {
      stop(
        "A selected full-network edge is ineligible in its fold: ", slug,
        " / ", replicate
      )
    }
  }
  full[, `:=`(
    timepoint = timepoint,
    eligible_folds = rowSums(eligible_matrix),
    selected_folds = rowSums(selected_matrix)
  )]
  full[, conditional_selection_fraction := fifelse(
    eligible_folds > 0L, selected_folds / eligible_folds, NA_real_
  )]
  full[, eligible_all_folds := eligible_folds == length(replicates)]
  full[, selected_all_eligible_folds :=
    eligible_folds > 0L & selected_folds == eligible_folds]

  summarize_subset <- function(x, top_k) {
    data.frame(
      dataset = "GSE274113", timepoint = timepoint, top_k = top_k,
      edges = nrow(x), folds = length(replicates),
      eligible_in_any_fold_fraction = mean(x$eligible_folds > 0L),
      eligible_in_all_folds_fraction = mean(x$eligible_all_folds),
      mean_conditional_selection_fraction = mean(
        x$conditional_selection_fraction,
        na.rm = TRUE
      ),
      median_conditional_selection_fraction = median(
        x$conditional_selection_fraction,
        na.rm = TRUE
      ),
      selected_in_all_eligible_folds_fraction = mean(
        x$selected_all_eligible_folds[x$eligible_folds > 0L]
      ),
      selected_in_all_folds_among_all_fold_eligible = mean(
        x$selected_folds[x$eligible_all_folds] == length(replicates)
      ), stringsAsFactors = FALSE
    )
  }
  top_rows[[slug]] <- rbindlist(lapply(
    unique(pmin(top_k_requested, nrow(full))),
    function(k) summarize_subset(full[seq_len(k)], k)
  ))
  decile_rows[[slug]] <- full[, .(
    edges = .N,
    eligible_in_any_fold_fraction = mean(eligible_folds > 0L),
    eligible_in_all_folds_fraction = mean(eligible_all_folds),
    mean_conditional_selection_fraction = mean(
      conditional_selection_fraction,
      na.rm = TRUE
    ),
    median_conditional_selection_fraction = median(
      conditional_selection_fraction,
      na.rm = TRUE
    ),
    selected_in_all_eligible_folds_fraction = mean(
      selected_all_eligible_folds[eligible_folds > 0L]
    ),
    selected_in_all_folds_among_all_fold_eligible = mean(
      selected_folds[eligible_all_folds] == length(replicates)
    )
  ), by = .(timepoint, rank_decile)]
  edge_rows[[slug]] <- full[, .(
    timepoint, regulator, target, weight, full_rank, rank_decile,
    target_detected_control_cells, regulator_detected_control_cells,
    eligible_folds, selected_folds, conditional_selection_fraction,
    eligible_all_folds, selected_all_eligible_folds
  )]
}

edges <- rbindlist(edge_rows)
top <- rbindlist(top_rows)
deciles <- rbindlist(decile_rows)
if (any(edges$selected_folds > edges$eligible_folds) ||
  any(!is.finite(top$mean_conditional_selection_fraction))) {
  stop("Eligibility-conditioned output contract failed.")
}
fwrite(edges, file.path(stage_root, "full_edge_conditioned_stability.tsv.gz"),
  sep = "\t"
)
fwrite(top, file.path(stage_root, "top_rank_conditioned_stability.tsv"), sep = "\t")
fwrite(deciles, file.path(stage_root, "rank_decile_conditioned_stability.tsv"),
  sep = "\t"
)
if (!file.rename(stage_root, final_root)) stop("Unable to promote conditioned output.")
cat("gse274113_eligibility_conditioned_stability\tPASS\t", final_root, "\n", sep = "")
