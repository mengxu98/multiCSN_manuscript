#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(multiCSN)
})

threads <- 4L  # Worker threads
if (!is.finite(threads) || threads < 1L) threads <- 4L
data.table::setDTthreads(threads)


# ---- Parameters ----
arguments <- list(
  "data/gse194122/gse194122_candidates.rds",
  "results/gse194122/trajectory/erythroid_coordinates_and_pseudotime.tsv.gz",
  "results/gse194122/trajectory"
)
# ---------------------------------------------------------------
input <- arguments[[1L]]
trajectory_path <- arguments[[2L]]
output <- arguments[[3L]]
input <- normalizePath(input, mustWork = TRUE)
trajectory_path <- normalizePath(trajectory_path, mustWork = TRUE)
stage <- tempfile(".result-files-", tmpdir = dirname(output))
if (file.exists(stage)) stop("Staging path already exists: ", stage, call. = FALSE)
dir.create(stage, recursive = TRUE)

trajectory <- fread(trajectory_path, select = c("cell_id", "state", "donor_id"))
if (!nrow(trajectory) || anyDuplicated(trajectory$cell_id) || anyNA(trajectory)) {
  stop("Invalid trajectory cell identifiers or metadata.")
}
states <- c("HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast")
trajectory[, state := factor(state, levels = states)]
setorder(trajectory, state, donor_id, cell_id)

object <- readRDS(input)
params <- multiCSN::Params(object)
gene_layer <- multiCSN::get_layer_data(object, assay = params$rna_assay, layer = "data")
requested_targets <- multiCSN::get_attribute(object, celltypes = "aggregate", attribute = "genes")
requested_tfs <- multiCSN::get_attribute(object, attribute = "tfs")
features <- sort(intersect(union(requested_targets, requested_tfs), rownames(gene_layer)), method = "radix")
if (!length(features)) stop("No expression features.")

groups <- unique(trajectory[, .(state, donor_id)])
groups[, group_id := seq_len(.N)]
trajectory <- merge(trajectory, groups, by = c("state", "donor_id"), sort = FALSE)
setorder(trajectory, group_id, cell_id)
cells <- trajectory$cell_id
group_design <- sparseMatrix(
  i = trajectory$group_id,
  j = seq_len(nrow(trajectory)),
  x = 1,
  dims = c(nrow(groups), nrow(trajectory))
)
group_cells <- as.integer(Matrix::rowSums(group_design))
expression <- gene_layer[features, cells, drop = FALSE]
group_sum <- expression %*% Matrix::t(group_design)
group_detected <- (expression != 0) %*% Matrix::t(group_design)
group_mean <- sweep(as.matrix(group_sum), 2L, group_cells, "/")
group_detection <- sweep(as.matrix(group_detected), 2L, group_cells, "/")

long_expression <- data.table(
  feature = rep(features, times = nrow(groups)),
  group_id = rep(seq_len(nrow(groups)), each = length(features)),
  mean_normalized_expression = as.vector(group_mean),
  detection_fraction = as.vector(group_detection)
)
long_expression <- merge(long_expression, groups, by = "group_id", sort = FALSE)
long_expression[, feature_role := fifelse(
  feature %in% requested_tfs & feature %in% requested_targets, "TF_and_target",
  fifelse(feature %in% requested_tfs, "TF", "target")
)]
setcolorder(long_expression, c(
  "state", "donor_id", "group_id", "feature", "feature_role",
  "mean_normalized_expression", "detection_fraction"
))
setorder(long_expression, state, donor_id, feature)

state_expression <- long_expression[, .(
  donor_groups = .N,
  median_donor_expression = median(mean_normalized_expression),
  q25_donor_expression = quantile(mean_normalized_expression, 0.25),
  q75_donor_expression = quantile(mean_normalized_expression, 0.75),
  median_detection_fraction = median(detection_fraction)
), by = .(state, feature, feature_role)]
state_expression[, state_order := match(state, states)]
setcolorder(state_expression, c("state_order", setdiff(names(state_expression), "state_order")))
setorder(state_expression, state_order, feature)

fwrite(long_expression, file.path(stage, "donor_state_expression.tsv.gz"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(state_expression, file.path(stage, "state_expression_summary.tsv.gz"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(groups[, .(
  state, donor_id, group_id,
  cells = group_cells[group_id]
)], file.path(stage, "donor_state_group_counts.tsv"), sep = "\t", quote = FALSE)
publish_result_files(stage, output)
