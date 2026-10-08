#!/usr/bin/env Rscript


suppressPackageStartupMessages(library(data.table))

# ---- Parameters ----
input_dir <- normalizePath("results/gse162170/loci/extracted", mustWork = TRUE)
output_dir <- "results/gse162170/loci/candidates"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(stem) {
  fread(file.path(input_dir, paste0(stem, ".tsv")), showProgress = FALSE)
}

triplets <- read_table("formal_focal_triplets")
accessibility <- read_table("focal_accessibility")
expression <- read_table("focal_expression")
expression_samples <- read_table("focal_expression_samples")
expression_cells <- read_table("focal_expression_cells")
accessibility_samples <- read_table("focal_accessibility_samples")
annotation <- read_table("focal_gene_annotation")
join_check <- read_table("exact_join_check")

stopifnot(
  nrow(join_check) == nrow(triplets),
  max(join_check$component_max_abs_delta, na.rm = TRUE) < 1e-12
)

parse_regions <- function(x) {
  pieces <- tstrsplit(sub("^chr", "", x$region), "-", fixed = TRUE)
  x[, `:=`(
    chromosome = paste0("chr", pieces[[1L]]),
    start = as.integer(pieces[[2L]]),
    end = as.integer(pieces[[3L]])
  )]
  x[]
}

targets <- sort(unique(triplets$target))
for (target_name in targets) {
  paths <- parse_regions(copy(triplets[target == target_name]))
  setorder(paths, start, end, -chain_delta_bic, regulator)

  acc <- copy(accessibility[region %chin% paths$region])
  acc <- parse_regions(acc)
  setorder(acc, start, end, broad_class)

  expr <- copy(expression[target == target_name])
  expr_samples <- copy(expression_samples[target == target_name])
  expr_cells <- copy(expression_cells[target == target_name])
  acc_samples <- copy(accessibility_samples[region %chin% paths$region])
  gene <- copy(annotation[gene_name == target_name])

  fwrite(paths, file.path(output_dir, paste0(target_name, "_paths.tsv")), sep = "\t")
  fwrite(acc, file.path(output_dir, paste0(target_name, "_accessibility.tsv")), sep = "\t")
  fwrite(acc_samples, file.path(output_dir, paste0(target_name, "_accessibility_samples.tsv")), sep = "\t")
  fwrite(expr, file.path(output_dir, paste0(target_name, "_expression_summary.tsv")), sep = "\t")
  fwrite(expr_samples, file.path(output_dir, paste0(target_name, "_expression_samples.tsv")), sep = "\t")
  fwrite(expr_cells, file.path(output_dir, paste0(target_name, "_expression_cells.tsv")), sep = "\t")
  fwrite(gene, file.path(output_dir, paste0(target_name, "_gene_annotation.tsv")), sep = "\t")
}

fwrite(
  triplets[, .(
    chains = .N,
    regions = uniqueN(region),
    regulators = uniqueN(regulator),
    positive = sum(chain_direction > 0),
    negative = sum(chain_direction < 0),
    median_chain_delta_bic = median(chain_delta_bic),
    q75_chain_delta_bic = quantile(chain_delta_bic, 0.75),
    max_chain_delta_bic = max(chain_delta_bic),
    backbone_selected = sum(backbone_selected),
    backbone_sign_concordant = sum(backbone_sign_concordant)
  ), by = target][order(target)],
  file.path(output_dir, "candidate_summary.tsv"),
  sep = "\t"
)
