#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(AUCell)
})

# ---- Parameters ----
object_path <- normalizePath("data/gse162170/gse162170_candidates.rds", mustWork = TRUE)
result_root <- normalizePath("results/gse162170/network", mustWork = TRUE)
out_dir <- normalizePath("results/gse162170/cell_programs", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
triplet_path <- file.path(result_root, "triplets.tsv")

group_mean <- function(mat, group, levels_fixed) {
  stopifnot(ncol(mat) == length(group), all(group %in% levels_fixed))
  design <- sparseMatrix(
    i = seq_along(group), j = match(group, levels_fixed), x = 1,
    dims = c(length(group), length(levels_fixed)),
    dimnames = list(colnames(mat), levels_fixed)
  )
  out <- as.matrix(mat %*% design)
  denom <- tabulate(match(group, levels_fixed), nbins = length(levels_fixed))
  stopifnot(all(denom > 0L))
  sweep(out, 2, denom, "/")
}

top_summary <- function(mat, feature_name) {
  stopifnot(nrow(mat) > 0L, ncol(mat) >= 2L)
  top_index <- max.col(mat, ties.method = "first")
  ordered <- t(apply(mat, 1, sort, decreasing = TRUE))
  out <- data.table(
    feature = rownames(mat),
    top_class = colnames(mat)[top_index],
    top_mean = mat[cbind(seq_len(nrow(mat)), top_index)],
    second_mean = ordered[, 2L]
  )
  out[, margin := top_mean - second_mean]
  setnames(out, "feature", feature_name)
  out
}

row_z <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x)) / s
}

message("Reading complete paired object and current exact-join triplets")
object <- readRDS(object_path)
triplets <- fread(triplet_path, showProgress = FALSE)

class_levels <- c("GluN", "IN", "IPC_RG", "OPC_Glia", "Vascular")
class_map <- c(
  "GluN2" = "GluN", "GluN3" = "GluN", "GluN4" = "GluN", "GluN5" = "GluN",
  "nIPC/GluN1" = "GluN", "SP" = "GluN",
  "IN1" = "IN", "IN2" = "IN", "IN3" = "IN",
  "RG" = "IPC_RG", "Cyc. Prog." = "IPC_RG",
  "mGPC/OPC" = "OPC_Glia", "EC/Peric." = "Vascular"
)

meta <- as.data.table(object[[]], keep.rownames = "cell")
meta[, broad_class := unname(class_map[author_cluster_name])]
stopifnot(
  nrow(meta) == 8981L,
  uniqueN(meta$cell) == 8981L,
  !anyNA(meta$broad_class),
  setequal(unique(meta$broad_class), class_levels),
  nrow(triplets) == 140786L,
  uniqueN(triplets$target) == 11017L,
  uniqueN(triplets$region) == 32969L,
  uniqueN(triplets$regulator) == 1038L,
  all(is.finite(triplets$chain_delta_bic)),
  all(triplets$chain_direction %in% c(-1, 1))
)

reduction_priority <- c("wnn.umap", "umap", "harmony.umap", "lsi")
available_reductions <- names(object@reductions)
selected_reduction <- reduction_priority[reduction_priority %in% available_reductions][1L]
if (is.na(selected_reduction)) {
  message("No stored 2D reduction; recomputing a visualization-only RNA UMAP from all cells")
  rna_for_embedding <- GetAssayData(object, assay = "RNA", layer = "data")
  n_cells_embedding <- ncol(rna_for_embedding)
  embedding_cells <- colnames(rna_for_embedding)
  stopifnot(length(embedding_cells) == 8981L, uniqueN(embedding_cells) == 8981L)
  gene_mean <- Matrix::rowMeans(rna_for_embedding)
  squared <- rna_for_embedding
  squared@x <- squared@x^2
  gene_variance <- (Matrix::rowSums(squared) - n_cells_embedding * gene_mean^2) /
    (n_cells_embedding - 1)
  rm(squared)
  gc(verbose = FALSE)
  variance_table <- data.table(
    gene = rownames(rna_for_embedding),
    variance = as.numeric(gene_variance)
  )
  variance_table <- variance_table[is.finite(variance) & variance > 0]
  setorder(variance_table, -variance, gene)
  variance_table[, variance_rank := seq_len(.N)]
  embedding_features <- variance_table[seq_len(min(2000L, .N))]
  stopifnot(nrow(embedding_features) == 2000L)
  write_tsv(embedding_features, file.path(out_dir, "GSE4a_embedding_features.tsv"))
  set.seed(20260815)
  pca <- irlba::prcomp_irlba(
    t(rna_for_embedding[embedding_features$gene, , drop = FALSE]),
    n = 30L, center = TRUE, scale. = FALSE
  )
  stopifnot(nrow(pca$x) == 8981L, ncol(pca$x) == 30L)
  set.seed(20260815)
  embedding <- uwot::umap(
    pca$x,
    n_neighbors = 30L,
    min_dist = 0.30,
    metric = "cosine",
    n_threads = 1L,
    n_sgd_threads = 1L,
    init = "spectral",
    verbose = FALSE
  )
  rownames(embedding) <- embedding_cells
  selected_reduction <- "recomputed_rna_umap_all_8981_cells"
  write_tsv(data.table(
    item = c(
      "stored_reduction_available", "visualization_recomputed", "cells",
      "variable_genes", "PCA_dimensions", "UMAP_neighbors",
      "UMAP_min_dist", "UMAP_metric", "seed", "cell_sampling",
      "Greedy_L0_refit"
    ),
    value = c(
      "false", "true", "8981", "2000", "30", "30", "0.30",
      "cosine", "20260815", "none", "none"
    )
  ), file.path(out_dir, "GSE4a_embedding_method.tsv"))
  rm(rna_for_embedding, pca)
  gc(verbose = FALSE)
} else {
  embedding <- Embeddings(object[[selected_reduction]])
  write_tsv(data.table(
    item = c("stored_reduction_available", "visualization_recomputed", "stored_reduction", "cells", "cell_sampling", "Greedy_L0_refit"),
    value = c("true", "false", selected_reduction, "8981", "none", "none")
  ), file.path(out_dir, "GSE4a_embedding_method.tsv"))
}
stopifnot(nrow(embedding) == nrow(meta), ncol(embedding) >= 2L)
embedding_dt <- data.table(
  cell = rownames(embedding),
  dimension_1 = embedding[, 1L],
  dimension_2 = embedding[, 2L]
)
embedding_dt <- merge(
  embedding_dt,
  meta[, .(cell, broad_class, author_cluster_name, Sample.ID, Sample.Age)],
  by = "cell", all.x = TRUE, sort = FALSE
)
stopifnot(nrow(embedding_dt) == 8981L, !anyNA(embedding_dt$broad_class))
composition <- meta[, .(
  cells = .N,
  fraction = .N / nrow(meta),
  samples = uniqueN(Sample.ID)
), by = .(broad_class, author_cluster_name)][order(match(broad_class, class_levels), -cells, author_cluster_name)]
class_composition <- meta[, .(
  cells = .N,
  fraction = .N / nrow(meta),
  samples = uniqueN(Sample.ID)
), by = broad_class][order(match(broad_class, class_levels))]
write_tsv(embedding_dt, file.path(out_dir, "GSE4a_cell_coordinates_and_labels.tsv"))
write_tsv(composition, file.path(out_dir, "GSE4a_cluster_composition.tsv"))
write_tsv(class_composition, file.path(out_dir, "GSE4a_class_composition.tsv"))

rna <- GetAssayData(object, assay = "RNA", layer = "data")
atac <- GetAssayData(object, assay = "ATAC", layer = "data")
stopifnot(ncol(rna) == 8981L, ncol(atac) == 8981L)

needed_genes <- intersect(unique(c(triplets$regulator, triplets$target)), rownames(rna))
needed_regions <- intersect(unique(triplets$region), rownames(atac))
rna_mean <- group_mean(rna[needed_genes, , drop = FALSE], meta$broad_class, class_levels)
atac_mean <- group_mean(atac[needed_regions, , drop = FALSE], meta$broad_class, class_levels)
rna_top <- top_summary(rna_mean, "gene")
atac_top <- top_summary(atac_mean, "region")

annotated <- merge(
  triplets,
  rna_top[, .(regulator = gene, regulator_top_class = top_class, regulator_margin = margin)],
  by = "regulator", all.x = TRUE
)
annotated <- merge(
  annotated,
  rna_top[, .(target = gene, target_top_class = top_class, target_margin = margin)],
  by = "target", all.x = TRUE
)
annotated <- merge(
  annotated,
  atac_top[, .(region, region_top_class = top_class, region_margin = margin)],
  by = "region", all.x = TRUE
)
stopifnot(!anyNA(annotated$target_top_class), !anyNA(annotated$region_top_class))
annotated[, `:=`(
  target_region_match = target_top_class == region_top_class,
  all_three_match = regulator_top_class == target_top_class & target_top_class == region_top_class,
  regulation_sign = fifelse(chain_direction >= 0, "positive", "negative")
)]

region_target <- unique(annotated[, .(
  target, region, region_top_class, region_margin, target_top_class, target_margin
)])
target_votes <- region_target[, .(selected_regions = .N), by = .(target, region_top_class)]
target_votes[, maximum := max(selected_regions), by = target]
target_votes[, tied_maximum := sum(selected_regions == maximum) > 1L, by = target]
target_anchor <- target_votes[
  selected_regions == maximum & !tied_maximum,
  .(target, chromatin_anchor = region_top_class, selected_regions)
]
target_anchor <- merge(
  target_anchor,
  unique(region_target[, .(target, target_top_class, target_margin)]),
  by = "target"
)
target_anchor[, rna_matches_chromatin_anchor := target_top_class == chromatin_anchor]
stopifnot(
  target_anchor[chromatin_anchor == "GluN", .N] == 2271L,
  target_anchor[chromatin_anchor == "IN", .N] == 1396L,
  target_anchor[chromatin_anchor == "IPC_RG", .N] == 1629L,
  target_anchor[chromatin_anchor == "OPC_Glia", .N] == 1308L,
  target_anchor[chromatin_anchor == "Vascular", .N] == 2444L
)
write_tsv(target_anchor[order(match(chromatin_anchor, class_levels), target)], file.path(out_dir, "target_membership.tsv"))

message("Computing full-cell AUCell rankings and five current target programs")
set.seed(20260815)
rankings <- AUCell_buildRankings(rna, nCores = 1, plotStats = FALSE, verbose = FALSE)
current_gene_sets <- split(target_anchor$target, target_anchor$chromatin_anchor)
current_gene_sets <- lapply(current_gene_sets, function(x) intersect(unique(x), rownames(rna)))
current_gene_sets <- current_gene_sets[class_levels]
stopifnot(all(lengths(current_gene_sets) > 0L))
auc_max_rank <- ceiling(0.05 * nrow(rna))
target_auc <- AUCell_calcAUC(
  current_gene_sets, rankings,
  nCores = 1,
  aucMaxRank = auc_max_rank, verbose = FALSE
)
target_auc_dt <- as.data.table(t(getAUC(target_auc)), keep.rownames = "cell")
target_auc_long <- melt(
  target_auc_dt,
  id.vars = "cell",
  variable.name = "target_program", value.name = "AUC"
)
target_auc_long <- merge(
  target_auc_long,
  meta[, .(cell, broad_class, author_cluster_name, Sample.ID)],
  by = "cell", all.x = TRUE
)
target_auc_summary <- target_auc_long[, .(
  mean_AUC = mean(AUC),
  median_AUC = median(AUC),
  cells = .N,
  samples = uniqueN(Sample.ID)
), by = .(target_program, broad_class)]
target_auc_summary[, top_class := broad_class[which.max(mean_AUC)], by = target_program]
target_auc_summary[, matching_class := target_program == broad_class]
stopifnot(
  all(unique(target_auc_summary[, .(target_program, top_class)])$target_program ==
    unique(target_auc_summary[, .(target_program, top_class)])$top_class),
  nrow(target_auc_long) == 8981L * 5L
)
write_tsv(target_auc_long[order(cell, match(target_program, class_levels))], file.path(out_dir, "aucell_by_cell.tsv.gz"))
write_tsv(target_auc_summary[order(match(target_program, class_levels), match(broad_class, class_levels))], file.path(out_dir, "aucell_by_class.tsv"))

module_summary <- annotated[, .(
  triplets = .N,
  targets = uniqueN(target),
  regions = uniqueN(region),
  target_region_match_fraction = mean(target_region_match),
  all_three_match_fraction = mean(all_three_match),
  median_chain_delta_bic = median(chain_delta_bic)
), by = .(regulator, anchor_class = region_top_class, sign = regulation_sign)]
setorder(module_summary, anchor_class, -targets, -regions, -triplets, regulator, sign)
module_summary[, rank := seq_len(.N), by = anchor_class]
top_modules <- module_summary[rank <= 5L]
top_modules[, module_id := paste(regulator, anchor_class, sign, sep = "__")]
stopifnot(nrow(top_modules) == 25L, all(top_modules[, .N, by = anchor_class]$N == 5L))

module_membership <- merge(
  annotated,
  top_modules[, .(regulator, anchor_class, sign, module_id, module_rank = rank)],
  by.x = c("regulator", "region_top_class", "regulation_sign"),
  by.y = c("regulator", "anchor_class", "sign"),
  all = FALSE
)
stopifnot(nrow(module_membership) > 0L)
membership_export <- module_membership[, .(
  module_id, module_rank, regulator,
  anchor_class = region_top_class,
  sign = regulation_sign, target, region,
  tf_region_delta_bic, region_gene_delta_bic,
  chain_delta_bic, chain_direction,
  target_top_class, region_top_class, regulator_top_class
)]
write_tsv(membership_export[order(match(anchor_class, class_levels), module_rank, target, region)], file.path(out_dir, "module_membership.tsv.gz"))
write_tsv(top_modules[order(match(anchor_class, class_levels), rank)], file.path(out_dir, "module_summary.tsv"))

message("Computing full-cell RNA AUCell and ATAC activity for 25 current modules")
module_targets <- split(module_membership$target, module_membership$module_id)
module_targets <- lapply(module_targets, function(x) intersect(unique(x), rownames(rna)))
module_regions <- split(module_membership$region, module_membership$module_id)
module_regions <- lapply(module_regions, function(x) intersect(unique(x), rownames(atac)))
module_ids <- top_modules[order(match(anchor_class, class_levels), rank), module_id]
module_targets <- module_targets[module_ids]
module_regions <- module_regions[module_ids]
stopifnot(all(lengths(module_targets) > 0L), all(lengths(module_regions) > 0L))

module_auc <- AUCell_calcAUC(
  module_targets, rankings,
  nCores = 1,
  aucMaxRank = auc_max_rank, verbose = FALSE
)
module_auc_mat <- t(getAUC(module_auc))
module_atac_mat <- vapply(module_regions, function(regions) {
  Matrix::colMeans(atac[regions, , drop = FALSE])
}, numeric(ncol(atac)))
rownames(module_atac_mat) <- colnames(atac)
stopifnot(
  nrow(module_auc_mat) == 8981L, ncol(module_auc_mat) == 25L,
  all(dim(module_atac_mat) == c(8981L, 25L)),
  setequal(rownames(module_auc_mat), meta$cell),
  setequal(rownames(module_atac_mat), meta$cell)
)

module_auc_dt <- as.data.table(module_auc_mat, keep.rownames = "cell")
module_auc_long <- melt(module_auc_dt, id.vars = "cell", variable.name = "module_id", value.name = "raw_activity")
module_auc_long[, modality := "RNA_target_AUCell"]
module_atac_dt <- as.data.table(module_atac_mat, keep.rownames = "cell")
module_atac_long <- melt(module_atac_dt, id.vars = "cell", variable.name = "module_id", value.name = "raw_activity")
module_atac_long[, modality := "ATAC_selected_region_mean"]
module_cell_activity <- rbindlist(list(module_auc_long, module_atac_long), use.names = TRUE)
module_cell_activity <- merge(
  module_cell_activity,
  meta[, .(cell, broad_class, author_cluster_name, Sample.ID)],
  by = "cell", all.x = TRUE
)
module_cell_activity <- merge(
  module_cell_activity,
  top_modules[, .(module_id, regulator, anchor_class, sign, module_rank = rank, targets, regions, triplets)],
  by = "module_id", all.x = TRUE
)
stopifnot(nrow(module_cell_activity) == 8981L * 25L * 2L, !anyNA(module_cell_activity$anchor_class))

module_raw_summary <- module_cell_activity[, .(
  mean_activity = mean(raw_activity),
  median_activity = median(raw_activity),
  cells = .N,
  samples = uniqueN(Sample.ID)
), by = .(module_id, regulator, anchor_class, sign, module_rank, targets, regions, triplets, modality, broad_class)]
module_raw_summary[, scaled_mean_within_module_modality := row_z(mean_activity), by = .(module_id, modality)]
module_raw_summary[, matching_class := anchor_class == broad_class]
write_tsv(module_cell_activity[order(module_id, modality, cell)], file.path(out_dir, "module_activity_by_cell.tsv.gz"))
write_tsv(module_raw_summary[order(match(anchor_class, class_levels), module_rank, modality, match(broad_class, class_levels))], file.path(out_dir, "module_raw_and_scaled_activity.tsv"))

concordance <- rbindlist(list(
  annotated[, .(
    unit = "triplet_rows", n = .N,
    target_region_match_n = sum(target_region_match),
    target_region_match_fraction = mean(target_region_match),
    all_three_match_n = sum(all_three_match),
    all_three_match_fraction = mean(all_three_match)
  )],
  unique(annotated[, .(target, region, target_region_match)])[, .(
    unit = "unique_region_target", n = .N,
    target_region_match_n = sum(target_region_match),
    target_region_match_fraction = mean(target_region_match),
    all_three_match_n = NA_integer_,
    all_three_match_fraction = NA_real_
  )]
), use.names = TRUE, fill = TRUE)
stopifnot(
  concordance[unit == "triplet_rows", n] == 140786L,
  concordance[unit == "triplet_rows", target_region_match_n] == 70983L,
  concordance[unit == "triplet_rows", all_three_match_n] == 32340L,
  concordance[unit == "unique_region_target", n] == 35158L,
  concordance[unit == "unique_region_target", target_region_match_n] == 16791L
)
write_tsv(concordance, file.path(out_dir, "supporting_cross_modality_concordance.tsv"))

message("Extraction complete: ", out_dir)
