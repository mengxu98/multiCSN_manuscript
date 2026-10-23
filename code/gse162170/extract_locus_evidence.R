#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(Signac)
})

# ---- Parameters ----
object_path <- normalizePath("data/gse162170/gse162170_candidates.rds", mustWork = TRUE)
formal_root <- normalizePath("results/gse162170/network", mustWork = TRUE)
output_dir <- "results/gse162170/loci/extracted"
locus_targets <- "ZNF804A,NFASC,KCNC2"
focal_targets <- strsplit(locus_targets, ",", fixed = TRUE)[[1L]]
focal_targets <- unique(trimws(focal_targets[nzchar(trimws(focal_targets))]))
if (!length(focal_targets)) {
  stop("TARGETS did not contain a target.", call. = FALSE)
}

emit <- function(name, x) {
  if (nzchar(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(x, file.path(output_dir, paste0(name, ".tsv")),
      sep = "\t",
      quote = FALSE, na = "NA"
    )
    return(invisible(NULL))
  }
  cat("@@BEGIN\t", name, "\n", sep = "")
  fwrite(x, file = "", sep = "\t", quote = FALSE, na = "NA")
  cat("@@END\t", name, "\n", sep = "")
}

group_summary <- function(mat, group, feature_name) {
  levels <- c("GluN", "IN", "IPC_RG", "OPC_Glia", "Vascular")
  design <- sparseMatrix(
    i = seq_along(group), j = match(group, levels), x = 1,
    dims = c(length(group), length(levels)),
    dimnames = list(colnames(mat), levels)
  )
  sums <- as.matrix(mat %*% design)
  counts <- tabulate(match(group, levels), nbins = length(levels))
  means <- sweep(sums, 2, counts, "/")
  detected <- as.matrix((mat > 0) %*% design)
  fractions <- sweep(detected, 2, counts, "/")
  mean_dt <- as.data.table(means, keep.rownames = feature_name)
  fraction_dt <- as.data.table(fractions, keep.rownames = feature_name)
  mean_long <- melt(mean_dt,
    id.vars = feature_name, variable.name = "broad_class",
    value.name = "mean_normalized_signal"
  )
  fraction_long <- melt(fraction_dt,
    id.vars = feature_name, variable.name = "broad_class",
    value.name = "detected_fraction"
  )
  merge(mean_long, fraction_long, by = c(feature_name, "broad_class"))
}

sample_summary <- function(mat, meta, feature_name) {
  groups <- unique(meta[, .(sample_id = Sample.ID, broad_class)])
  setorder(groups, sample_id, broad_class)
  group_id <- match(
    paste(meta$Sample.ID, meta$broad_class, sep = "\r"),
    paste(groups$sample_id, groups$broad_class, sep = "\r")
  )
  design <- sparseMatrix(
    i = seq_len(nrow(meta)), j = group_id, x = 1,
    dims = c(nrow(meta), nrow(groups)),
    dimnames = list(meta$cell, NULL)
  )
  counts <- tabulate(group_id, nbins = nrow(groups))
  means <- sweep(as.matrix(mat %*% design), 2, counts, "/")
  detected <- sweep(as.matrix((mat > 0) %*% design), 2, counts, "/")
  out <- rbindlist(lapply(seq_len(nrow(groups)), function(i) {
    data.table(
      feature = rownames(mat),
      sample_id = groups$sample_id[[i]],
      broad_class = groups$broad_class[[i]],
      cells = counts[[i]],
      mean_normalized_signal = means[, i],
      detected_fraction = detected[, i]
    )
  }))
  setnames(out, "feature", feature_name)
  out
}

triplets <- fread(file.path(formal_root, "triplets.tsv"), showProgress = FALSE)
region_gene <- fread(file.path(formal_root, "region_gene.tsv"), showProgress = FALSE)
tf_region <- fread(file.path(formal_root, "tf_region.tsv"), showProgress = FALSE)
focal_triplets <- triplets[target %chin% focal_targets]
focal_region_gene <- region_gene[target %chin% focal_targets]
focal_regions <- unique(focal_triplets$region)
focal_tf_region <- tf_region[region %chin% focal_regions]

setnames(
  focal_tf_region,
  c("standardized_beta", "deletion_delta_bic", "weight"),
  c("tf_region_beta", "tf_region_delta_bic", "tf_region_weight")
)
setnames(
  focal_region_gene,
  c("standardized_beta", "deletion_delta_bic", "weight"),
  c("region_gene_beta", "region_gene_delta_bic", "region_gene_weight")
)
exact_join <- merge(
  focal_tf_region, focal_region_gene,
  by = "region", allow.cartesian = TRUE, sort = FALSE
)[target %chin% focal_targets]
exact_join[, `:=`(
  chain_delta_bic_recomputed = pmin(tf_region_delta_bic, region_gene_delta_bic),
  chain_direction_recomputed = sign(tf_region_beta * region_gene_beta)
)]
setcolorder(exact_join, c(
  "regulator", "target", "region", "tf_region_beta", "tf_region_delta_bic",
  "tf_region_weight", "region_gene_beta", "region_gene_delta_bic",
  "region_gene_weight", "chain_delta_bic_recomputed", "chain_direction_recomputed"
))

join_check <- merge(
  focal_triplets,
  exact_join,
  by = c("regulator", "target", "region"),
  all = TRUE, suffixes = c("_formal_triplet", "_exact_join")
)
join_check[, component_max_abs_delta := pmax(
  abs(tf_region_beta_formal_triplet - tf_region_beta_exact_join),
  abs(tf_region_delta_bic_formal_triplet - tf_region_delta_bic_exact_join),
  abs(region_gene_beta_formal_triplet - region_gene_beta_exact_join),
  abs(region_gene_delta_bic_formal_triplet - region_gene_delta_bic_exact_join),
  abs(chain_delta_bic - chain_delta_bic_recomputed),
  abs(chain_direction - chain_direction_recomputed),
  na.rm = TRUE
)]

object <- readRDS(object_path)
meta <- as.data.table(object[[]], keep.rownames = "cell")
class_map <- c(
  "GluN2" = "GluN", "GluN3" = "GluN", "GluN4" = "GluN", "GluN5" = "GluN",
  "nIPC/GluN1" = "GluN", "SP" = "GluN",
  "IN1" = "IN", "IN2" = "IN", "IN3" = "IN",
  "RG" = "IPC_RG", "Cyc. Prog." = "IPC_RG",
  "mGPC/OPC" = "OPC_Glia", "EC/Peric." = "Vascular"
)
meta[, broad_class := unname(class_map[author_cluster_name])]

rna <- GetAssayData(object, assay = "RNA", layer = "data")
atac <- GetAssayData(object, assay = "ATAC", layer = "data")
expression <- group_summary(rna[focal_targets, , drop = FALSE], meta$broad_class, "target")
accessibility <- group_summary(atac[focal_regions, , drop = FALSE], meta$broad_class, "region")
expression_samples <- sample_summary(
  rna[focal_targets, , drop = FALSE], meta, "target"
)
accessibility_samples <- sample_summary(
  atac[focal_regions, , drop = FALSE], meta, "region"
)
expression_cells <- as.data.table(
  as.matrix(t(rna[focal_targets, , drop = FALSE])),
  keep.rownames = "cell"
)
expression_cells <- merge(
  meta[, .(cell, author_cluster_name, broad_class, sample_id = Sample.ID)],
  expression_cells,
  by = "cell", all.x = TRUE, sort = FALSE
)
expression_cells <- melt(
  expression_cells,
  id.vars = c("cell", "author_cluster_name", "broad_class", "sample_id"),
  variable.name = "target", value.name = "normalized_expression"
)

annotation <- as.data.table(as.data.frame(Signac::Annotation(object[["ATAC"]])))
annotation[, seqnames := as.character(seqnames)]
annotation[, strand := as.character(strand)]
gene_column <- intersect(c("gene_name", "gene", "symbol"), names(annotation))[[1L]]
focal_annotation <- annotation[get(gene_column) %chin% focal_targets]
if (!nrow(focal_annotation)) stop("Focal genes are absent from the ATAC annotation.")
annotation_ranges <- rbindlist(lapply(focal_targets, function(target_name) {
  target_rows <- focal_annotation[get(gene_column) == target_name]
  chromosome <- unique(target_rows$seqnames)
  target_regions <- focal_triplets[target == target_name, region]
  parsed <- tstrsplit(sub("^chr", "", target_regions), "-", fixed = TRUE)
  region_chromosome <- paste0("chr", parsed[[1L]])
  region_start <- as.integer(parsed[[2L]])
  region_end <- as.integer(parsed[[3L]])
  lower <- min(c(target_rows$start, region_start), na.rm = TRUE) - 10000L
  upper <- max(c(target_rows$end, region_end), na.rm = TRUE) + 10000L
  annotation[seqnames %chin% unique(c(chromosome, region_chromosome)) &
    end >= lower & start <= upper][, focal_target := target_name]
}), use.names = TRUE, fill = TRUE)

composition <- meta[, .(
  cells = .N,
  samples = uniqueN(Sample.ID)
), by = .(broad_class)][match(c("GluN", "IN", "IPC_RG", "OPC_Glia", "Vascular"), broad_class)]

emit("formal_focal_triplets", focal_triplets)
emit("formal_focal_region_gene", focal_region_gene)
emit("formal_focal_tf_region", focal_tf_region)
emit("exact_component_join", exact_join)
emit("exact_join_check", join_check)
emit("class_composition", composition)
emit("focal_expression", expression)
emit("focal_expression_samples", expression_samples)
emit("focal_expression_cells", expression_cells)
emit("focal_accessibility", accessibility)
emit("focal_accessibility_samples", accessibility_samples)
emit("focal_gene_annotation", focal_annotation)
emit("locus_annotation", annotation_ranges)
cat("@@SESSION\n")
cat("@@COMPLETE\n")
