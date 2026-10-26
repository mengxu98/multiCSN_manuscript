#!/usr/bin/env Rscript

source("code/functions/network_io.R")

source("code/functions/runtime.R")
configure_task_runtime()

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rhdf5)
  library(edgeR)
  library(data.table)
  library(readxl)
})
data.table::setDTthreads(1L)

# ---- Parameters ----
script_options <- list(cores = 1L, output = NULL)
cores <- script_options$cores
# ---------------------------------------------------------------
if (!identical(cores, 1L)) stop("This benchmark requires cores = 1.")
repo_root <- normalizePath(getwd(), mustWork = TRUE)
data_root <- file.path(repo_root, "data", "gse274113", "GSE274113")
inventory_root <- file.path(
  repo_root, "results", "gse274113", "inventory"
)
destination <- if (is.null(script_options$output)) {
  file.path(repo_root, "data", "gse274113", "heldout_pseudobulk_full")
} else {
  normalizePath(script_options$output, mustWork = FALSE)
}
output_parent <- dirname(destination)
timepoints <- c("day 7", "day 9", "day 11", "day 14")
perturbation_type <- "coding-sequence loss-of-function"
expected_function_direction <- "reduced TF function"
direction_validation_basis <- paste(
  "study-reported coding-sequence loss-of-function design and released target annotation;",
  "validates the intervention-direction contract, not per-guide perturbation efficacy"
)

metadata_path <- file.path(data_root, "GSE274113_annotated_metadata.csv.gz")
object_root <- file.path(
  repo_root, "data", "gse274113", "objects", "GSE274113_control_full"
)
object_counts_path <- file.path(object_root, "GSE274113_all_NT_controls_counts.h5")
object_validation_path <- file.path(object_root, "control_object_validation.tsv")
supplement_root <- file.path(data_root, "supplementary")
table_s2_path <- file.path(supplement_root, "NIHMS2076308-supplement-Table_S2.xlsx")
table_s3_path <- file.path(supplement_root, "NIHMS2076308-supplement-Table_S3.xlsx")
inventory_validation <- file.path(inventory_root, "validation.tsv")
if (dir.exists(destination)) stop("Formal causal truth already exists: ", destination)

inventory_values <- data.table::fread(inventory_validation)
if (!identical(colnames(inventory_values), c("field", "value")) ||
  anyDuplicated(inventory_values$field)) {
  stop("GSE274113 inventory validation schema is invalid.")
}
inventory_value <- function(field) {
  value <- inventory_values$value[inventory_values$field == field]
  if (length(value) != 1L) stop("Missing inventory field: ", field)
  value[[1]]
}

raw_metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
cell_column_index <- which(
  is.na(colnames(raw_metadata)) |
    !nzchar(colnames(raw_metadata)) |
    grepl("^Unnamed:", colnames(raw_metadata))
)
if (length(cell_column_index) != 1L) {
  stop("Unable to identify the unique released cell-ID column.")
}
released_cell_ids <- as.character(raw_metadata[[cell_column_index]])
informative <- !is.na(released_cell_ids) & nzchar(trimws(released_cell_ids))
metadata <- raw_metadata[informative, , drop = FALSE]
metadata$cell_id <- as.character(metadata[[cell_column_index]])
if (!nrow(metadata) || anyDuplicated(metadata$cell_id) ||
  !all(c("target", "replicate", "Timepoint") %in% names(metadata)) ||
  anyNA(metadata$target) || any(!nzchar(metadata$target)) ||
  !"replicate" %in% names(metadata) || anyNA(metadata$replicate) ||
  any(!grepl("^rep[1-9][0-9]*$", metadata$replicate))) {
  stop("Released metadata failed full-cell validation.")
}
targets <- sort(setdiff(unique(metadata$target), "NT"))
if (!length(targets) || !any(metadata$target == "NT")) {
  stop("Released perturbation/control split failed validation.")
}

replicates <- sort(unique(as.integer(sub("^rep", "", metadata$replicate))))
rna_guide_paths <- file.path(
  data_root, paste0("GSE274113_filtered_feature_bc_matrix_", replicates, ".h5")
)
multiome_paths <- file.path(
  data_root, paste0("GSE274113_rep", replicates, "_filtered_feature_bc_matrix.h5")
)

object_validation <- data.table::fread(object_validation_path)
if (!identical(colnames(object_validation), c("field", "value")) ||
  anyDuplicated(object_validation$field)) {
  stop("Control-object validation schema is invalid.")
}
expected_consensus_features <- as.integer(
  object_validation[field == "atac_consensus_features", value]
)
if (length(expected_consensus_features) != 1L || !is.finite(expected_consensus_features)) {
  stop("Control-object validation lacks a unique consensus-feature count.")
}
consensus_peak_names <- as.character(rhdf5::h5read(object_counts_path, "atac/features"))
if (length(consensus_peak_names) != expected_consensus_features ||
  anyDuplicated(consensus_peak_names)) {
  stop("Control-network consensus peak universe failed validation.")
}
consensus_peaks <- Signac::StringToGRanges(
  regions = consensus_peak_names,
  sep = c(":", "-")
)
standard_seqlevels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
stage <- tempfile(pattern = ".staging_GSE274113_causal_truth_", tmpdir = output_parent)
dir.create(stage)
dir.create(file.path(stage, "rna_effects"))
dir.create(file.path(stage, "atac_effects"))

rna_aggregates <- vector("list", length(replicates))
atac_aggregates <- vector("list", length(replicates))
sample_rows <- vector("list", length(replicates))
reference_rna_features <- NULL
nonstandard_peak_rows <- integer(length(replicates))
nonstandard_count_totals <- numeric(length(replicates))

for (index in seq_along(replicates)) {
  rep <- replicates[[index]]
  rep_name <- paste0("rep", rep)
  cat("aggregate_started\t", rep_name, "\t", format(Sys.time(), usetz = TRUE), "\n", sep = "")
  counts <- Seurat::Read10X_h5(
    rna_guide_paths[[index]],
    use.names = TRUE, unique.features = TRUE
  )
  if (!is.list(counts) ||
    !all(c("Gene Expression", "CRISPR Guide Capture") %in% names(counts))) {
    stop("Missing RNA or guide feature type for ", rep_name)
  }
  rna <- counts[["Gene Expression"]]
  if (is.null(reference_rna_features)) {
    reference_rna_features <- rownames(rna)
  } else if (!identical(rownames(rna), reference_rna_features)) {
    stop("Released RNA feature order differs across replicates at ", rep_name)
  }
  rep_metadata <- metadata[metadata$replicate == rep_name, , drop = FALSE]
  rep_metadata$raw_barcode <- unname(sub(
    paste0("^", rep_name, "_"), "", rep_metadata$cell_id
  ))
  if (anyDuplicated(rep_metadata$raw_barcode)) stop("Duplicate cell barcodes: ", rep_name)
  rna <- rna[, rep_metadata$raw_barcode, drop = FALSE]
  if (!all(as.numeric(Matrix::colSums(rna)) == as.numeric(rep_metadata$nCount_RNA))) {
    stop("Released RNA counts disagree with metadata for ", rep_name)
  }
  multiome_counts <- Seurat::Read10X_h5(
    multiome_paths[[index]],
    use.names = TRUE, unique.features = TRUE
  )
  if (!is.list(multiome_counts) ||
    !all(c("Gene Expression", "Peaks") %in% names(multiome_counts))) {
    stop("Released RNA+ATAC H5 lacks a peak matrix for ", rep_name)
  }
  release_atac <- multiome_counts[["Peaks"]]
  release_atac <- release_atac[, colnames(rna), drop = FALSE]
  if (!"nCount_atac" %in% colnames(rep_metadata)) {
    stop("Released metadata lacks nCount_atac for ", rep_name)
  }
  if (!all(as.numeric(Matrix::colSums(release_atac)) ==
    as.numeric(rep_metadata$nCount_atac))) {
    stop("Released ATAC counts disagree with metadata for ", rep_name)
  }
  release_ranges <- Signac::StringToGRanges(
    regions = rownames(release_atac), sep = c(":", "-")
  )
  standard_peaks <- as.character(GenomeInfoDb::seqnames(release_ranges)) %in%
    standard_seqlevels
  nonstandard_peak_rows[[index]] <- sum(!standard_peaks)
  nonstandard_count_totals[[index]] <- sum(
    release_atac[!standard_peaks, , drop = FALSE]
  )
  release_ranges <- release_ranges[standard_peaks]
  release_atac <- release_atac[standard_peaks, , drop = FALSE]
  within_replicate_overlaps <- GenomicRanges::findOverlaps(
    release_ranges, release_ranges,
    ignore.strand = TRUE
  )
  if (any(tabulate(
    S4Vectors::queryHits(within_replicate_overlaps),
    nbins = length(release_ranges)
  ) != 1L)) {
    stop(
      "Released peaks overlap within ", rep_name,
      "; consensus aggregation would double-count peak-matrix entries."
    )
  }
  # Consensus projection lives in the package (one-to-one mapping is validated there).
  consensus_index <- multiCSN::map_peaks_to_consensus(
    release_ranges, consensus_peaks, type = "within", require_one_to_one = TRUE
  )
  atac <- multiCSN::project_counts_to_consensus(
    release_atac, rownames(release_atac), consensus_peaks, type = "within"
  )
  if (!identical(rownames(atac), consensus_peak_names) ||
    !identical(colnames(atac), colnames(rna))) {
    stop("Consensus H5 peak aggregation changed features or cells for ", rep_name)
  }
  rep_timepoint <- unique(rep_metadata$Timepoint)
  if (length(rep_timepoint) != 1L || !rep_timepoint %in% timepoints) {
    stop("Replicate does not map to exactly one declared timepoint: ", rep_name)
  }
  # Paired replicate x condition aggregation lives in the package.
  aggregated <- multiCSN::pseudobulk_perturbation(
    list(RNA = rna, ATAC = atac),
    sample_info = data.frame(
      cell_id = colnames(rna),
      replicate = rep_name,
      timepoint = rep_timepoint,
      target = as.character(rep_metadata$target),
      stringsAsFactors = FALSE
    ),
    targets = targets, control_label = "NT", verbose = FALSE
  )
  if (sum(aggregated$sample_metadata$cells) != ncol(rna)) {
    stop("At least one perturbation group is absent from ", rep_name)
  }
  sample_rows[[index]] <- aggregated$sample_metadata
  rna_aggregates[[index]] <- aggregated$pseudobulk$RNA
  atac_aggregates[[index]] <- aggregated$pseudobulk$ATAC
  rm(
    counts, rna, atac, multiome_counts, release_atac, release_ranges,
    within_replicate_overlaps, consensus_index, aggregated
  )
  gc(verbose = FALSE)
  cat("aggregate_completed\t", rep_name, "\t", format(Sys.time(), usetz = TRUE), "\n", sep = "")
}

rna_counts <- do.call(cbind, rna_aggregates)
atac_counts <- do.call(cbind, atac_aggregates)
sample_metadata <- do.call(rbind, sample_rows)
rownames(sample_metadata) <- sample_metadata$sample_id
if (!identical(colnames(rna_counts), sample_metadata$sample_id) ||
  !identical(colnames(atac_counts), sample_metadata$sample_id) ||
  sum(sample_metadata$cells) != nrow(metadata) ||
  sum(sample_metadata$cells[sample_metadata$role == "network_training_control"]) != sum(metadata$target == "NT") ||
  sum(sample_metadata$cells[sample_metadata$role == "heldout_perturbation"]) != sum(metadata$target != "NT")) {
  stop("Pseudobulk aggregation failed complete-cell accounting.")
}
write_tsv(sample_metadata, file.path(stage, "pseudobulk_sample_metadata.tsv"), na = "")
saveRDS(rna_counts, file.path(stage, "rna_pseudobulk_counts.rds"), compress = FALSE)
saveRDS(atac_counts, file.path(stage, "atac_pseudobulk_counts.rds"), compress = FALSE)

# Paired perturbation models (edgeR quasi-likelihood with a replicate block)
# live in the package; the dataset contract below only relabels the summary.
effects <- multiCSN::fit_perturbation_effects(
  pseudobulk = list(RNA = rna_counts, ATAC = atac_counts),
  sample_metadata = sample_metadata,
  targets = targets, timepoints = timepoints,
  control_label = "NT",
  effect_definition = "heldout_perturbed_minus_NT_control",
  output_dir = stage,
  output_subdirs = c(RNA = "rna_effects", ATAC = "atac_effects"),
  verbose = TRUE
)
effect_summary <- effects$summary
effect_summary$output_file <- sub(paste0("^", stage, "/?"), "", effect_summary$output_file)
effect_summary$perturbation_type <- perturbation_type
effect_summary$expected_function_direction <- expected_function_direction
effect_summary$direction_validated <- TRUE
effect_summary$validation_basis <- direction_validation_basis
effect_summary <- effect_summary[, c(
  "modality", "timepoint", "target", "perturbation_type",
  "expected_function_direction", "direction_validated", "validation_basis",
  "paired_replicates", "control_cells", "perturbation_cells", "release_features",
  "tested_features", "significant_fdr_0_05", "significant_up_fdr_0_05",
  "significant_down_fdr_0_05", "output_file"
)]
expected_effect_design <- data.table::CJ(
  modality = c("ATAC", "RNA"),
  timepoint = timepoints,
  target = targets,
  unique = TRUE
)
observed_effect_design <- unique(data.table::as.data.table(effect_summary)[
  , .(modality, timepoint, target)
])
if (anyDuplicated(data.table::as.data.table(effect_summary)[, .(modality, timepoint, target)]) ||
  !data.table::fsetequal(observed_effect_design, expected_effect_design) ||
  any(effect_summary$perturbation_type != perturbation_type) ||
  any(effect_summary$expected_function_direction != expected_function_direction) ||
  any(!effect_summary$direction_validated) ||
  any(effect_summary$validation_basis != direction_validation_basis)) {
  stop("Perturbation-direction metadata or the input TF/timepoint truth design is incomplete.")
}
write_tsv(effect_summary, file.path(stage, "effect_model_summary.tsv"), na = "")

table_s2 <- as.data.frame(readxl::read_excel(table_s2_path, skip = 2))
table_s3 <- as.data.frame(readxl::read_excel(table_s3_path, skip = 2))
expected_columns <- c("log2FC", "beta_weight", "p_weight", "perturbation_name")
if (!identical(colnames(table_s2), c("peak_locus", expected_columns)) ||
  !identical(colnames(table_s3), c("gene_ID", expected_columns)) ||
  !nrow(table_s2) || !nrow(table_s3) ||
  anyNA(table_s2$peak_locus) || anyNA(table_s3$gene_ID)) {
  stop("Author S2/S3 tables differ from the audited release.")
}
write_tsv(table_s2, file.path(stage, "author_Table_S2_TF_sensitive_ACRs.tsv.gz"), na = "")
write_tsv(table_s3, file.path(stage, "author_Table_S3_TF_sensitive_genes.tsv.gz"), na = "")
author_summary <- data.frame(
  table = c("S2", "S3"),
  endpoint = c("TF-sensitive ACR", "TF-sensitive gene symbol"),
  rows = c(nrow(table_s2), nrow(table_s3)),
  unique_entities = c(length(unique(table_s2$peak_locus)), length(unique(table_s3$gene_ID))),
  perturbations = c(length(unique(table_s2$perturbation_name)), length(unique(table_s3$perturbation_name))),
  contains_timepoint_field = FALSE,
  contains_replicate_field = FALSE,
  permitted_evaluation_scope = "TF-level aggregate concordance only",
  stringsAsFactors = FALSE
)
write_tsv(author_summary, file.path(stage, "author_truth_scope.tsv"), na = "")

if (!file.rename(stage, destination)) stop("Unable to promote formal causal truth.")
print(effect_summary)
print(author_summary)
