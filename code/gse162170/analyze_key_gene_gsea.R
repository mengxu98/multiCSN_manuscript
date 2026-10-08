#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(fgsea)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
})
data.table::setDTthreads(1L)

# ---- Parameters ----
object_path <- normalizePath("data/gse162170/gse162170_candidates.rds", mustWork = TRUE)
result_root <- normalizePath("results/gse162170/network", mustWork = TRUE)
out_dir <- normalizePath("results/gse162170/key_gene_gsea", mustWork = FALSE)
# ---------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
triplet_path <- file.path(result_root, "triplets.tsv")
stopifnot(file.exists(triplet_path))

class_levels <- c("GluN", "IN", "IPC_RG", "OPC_Glia", "Vascular")
class_map <- c(
  "GluN2" = "GluN", "GluN3" = "GluN", "GluN4" = "GluN",
  "GluN5" = "GluN", "nIPC/GluN1" = "GluN", "SP" = "GluN",
  "IN1" = "IN", "IN2" = "IN", "IN3" = "IN",
  "RG" = "IPC_RG", "Cyc. Prog." = "IPC_RG",
  "mGPC/OPC" = "OPC_Glia", "EC/Peric." = "Vascular"
)

message("Reading the complete paired-cell object and exact-join targets")
object <- readRDS(object_path)
triplets <- fread(triplet_path, select = "target", showProgress = FALSE)
targets <- sort(unique(triplets$target))
meta <- as.data.table(object[[]], keep.rownames = "cell")
meta[, broad_class := unname(class_map[author_cluster_name])]
stopifnot(
  nrow(meta) == 8981L, uniqueN(meta$cell) == 8981L,
  !anyNA(meta$broad_class), setequal(unique(meta$broad_class), class_levels),
  uniqueN(meta$Sample.ID) == 3L, length(targets) == 11017L
)

rna <- GetAssayData(object, assay = "RNA", layer = "data")
targets <- intersect(targets, rownames(rna))
stopifnot(length(targets) == 11017L, identical(colnames(rna), meta$cell))
rna <- rna[targets, , drop = FALSE]

samples <- sort(unique(meta$Sample.ID))
group_levels <- as.vector(outer(samples, class_levels, paste, sep = "__"))
group <- paste(meta$Sample.ID, meta$broad_class, sep = "__")
design <- sparseMatrix(
  i = seq_along(group), j = match(group, group_levels), x = 1,
  dims = c(length(group), length(group_levels)),
  dimnames = list(meta$cell, group_levels)
)
group_sum <- as.matrix(rna %*% design)
group_n <- tabulate(match(group, group_levels), nbins = length(group_levels))
stopifnot(all(group_n > 0L))
group_mean <- sweep(group_sum, 2L, group_n, "/")

rank_rows <- lapply(class_levels, function(focal) {
  sample_effects <- vapply(samples, function(sample_id) {
    focal_group <- paste(sample_id, focal, sep = "__")
    other_groups <- paste(sample_id, setdiff(class_levels, focal), sep = "__")
    focal_mean <- group_mean[, focal_group]
    other_n <- group_n[match(other_groups, group_levels)]
    other_mean <- rowSums(group_sum[, other_groups, drop = FALSE]) / sum(other_n)
    focal_mean - other_mean
  }, numeric(length(targets)))
  colnames(sample_effects) <- samples
  score <- rowMeans(sample_effects)
  x <- data.table(
    anchor_class = focal,
    gene = targets,
    rank_score = score,
    samples_positive = rowSums(sample_effects > 0),
    samples_negative = rowSums(sample_effects < 0),
    min_sample_effect = apply(sample_effects, 1L, min),
    max_sample_effect = apply(sample_effects, 1L, max)
  )
  setorder(x, -rank_score, gene)
  x[, rank := seq_len(.N)]
  x
})
gene_ranks <- rbindlist(rank_rows)
stopifnot(
  nrow(gene_ranks) == length(class_levels) * length(targets),
  all(gene_ranks[, .N, by = anchor_class]$N == length(targets))
)
write_tsv(gene_ranks, file.path(out_dir, "gene_rank_scores.tsv.gz"))

message("Building GO biological-process pathways from the selected-target universe")
go_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = targets, keytype = "SYMBOL",
  columns = c("GOALL", "ONTOLOGYALL")
)
go_map <- as.data.table(go_map)[
  ONTOLOGYALL == "BP" & !is.na(GOALL) & SYMBOL %chin% targets,
  .(term_id = GOALL, gene = SYMBOL)
]
go_map <- unique(go_map)
pathways <- split(go_map$gene, go_map$term_id)
pathways <- pathways[lengths(pathways) >= 15L & lengths(pathways) <= 500L]
term_names <- AnnotationDbi::Term(GOTERM[names(pathways)])
stopifnot(length(pathways) > 1000L, !anyNA(term_names))

message("Running five sample-blocked preranked GSEA analyses")
set.seed(20260825)
gsea_rows <- lapply(class_levels, function(focal) {
  x <- copy(gene_ranks[anchor_class == focal])
  setorder(x, -rank_score, gene)
  tie_offset <- rev(seq_len(nrow(x))) * 1e-12
  stats <- setNames(x$rank_score + tie_offset, x$gene)
  fg <- fgseaMultilevel(
    pathways = pathways, stats = stats, minSize = 15L, maxSize = 500L,
    eps = 0, scoreType = "std", nPermSimple = 10000L
  )
  fg <- as.data.table(fg)
  fg[, `:=`(
    anchor_class = focal,
    Description = unname(term_names[pathway]),
    leading_edge_n = lengths(leadingEdge),
    leading_edge_genes = vapply(leadingEdge, paste, collapse = ";", character(1))
  )]
  fg[, leadingEdge := NULL]
  fg
})
gsea_all <- rbindlist(gsea_rows, use.names = TRUE, fill = TRUE)
gsea_all[, abs_NES_for_order := abs(NES)]
setorder(gsea_all, anchor_class, padj, pval, -abs_NES_for_order, pathway)
gsea_all[, abs_NES_for_order := NULL]
gsea_all[, class_rank := seq_len(.N), by = anchor_class]
gsea_all[, significant := is.finite(padj) & padj < 0.05]
write_tsv(gsea_all, file.path(out_dir, "gsea_all_terms.tsv.gz"))

gsea_display <- gsea_all[significant == TRUE & is.finite(NES)]
gsea_display[, abs_NES_for_order := abs(NES)]
setorder(gsea_display, anchor_class, padj, -abs_NES_for_order, pathway)
gsea_display <- gsea_display[, head(.SD, 3L), by = anchor_class]
gsea_display[, abs_NES_for_order := NULL]
gsea_display[, display_rank := seq_len(.N), by = anchor_class]
write_tsv(gsea_display, file.path(out_dir, "gsea_display_terms.tsv"))

gsea_display[, display_row_id := seq_len(.N)]
display_edges <- gsea_display[NES > 0, .(
  anchor_class, pathway, Description, NES, padj, significant,
  gene = unlist(strsplit(leading_edge_genes, ";", fixed = TRUE))
), by = display_row_id][, display_row_id := NULL]
display_edges <- merge(
  display_edges,
  gene_ranks[, .(anchor_class, gene, rank_score, rank, samples_positive)],
  by = c("anchor_class", "gene"), all.x = TRUE
)
key_genes <- display_edges[, .(
  enriched_terms = uniqueN(pathway),
  strongest_NES = max(NES),
  best_padj = min(padj, na.rm = TRUE),
  rank_score = unique(rank_score),
  rank = unique(rank),
  samples_positive = unique(samples_positive)
), by = .(anchor_class, gene)]
setorder(key_genes, anchor_class, -enriched_terms, -rank_score, gene)
key_genes[, key_gene_rank := seq_len(.N), by = anchor_class]
write_tsv(key_genes, file.path(out_dir, "gsea_key_genes.tsv"))

message("GSEA complete: ", out_dir)
