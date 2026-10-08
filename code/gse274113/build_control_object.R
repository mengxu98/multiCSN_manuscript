#!/usr/bin/env Rscript

source("code/functions/network_io.R")

source("code/functions/runtime.R")
configure_task_runtime()

replicates <- c(1:10, 12:14, 16)
repo_root <- normalizePath(getwd(), mustWork = TRUE)
script_path <- normalizePath("code/gse274113/build_control_object.R", mustWork = TRUE)
data_root <- file.path(repo_root, "data", "gse274113", "GSE274113")
extracted_root <- file.path(data_root, "extracted")
inventory_root <- file.path(
  repo_root, "results", "gse274113", "inventory", "GSE274113_full"
)
script_options <- list(output = NULL, cores = 1L)
if (!identical(script_options$cores, 1L)) stop("This build requires cores = 1.")
destination <- if (is.null(script_options$output)) {
  file.path(repo_root, "data", "gse274113", "objects", "GSE274113_control_full")
} else {
  normalizePath(script_options$output, mustWork = FALSE)
}
output_parent <- dirname(destination)
if (is.null(script_options$output)) {
  paths <- file.path(repo_root, "results/gse274113/networks",
    c("multiCSN", "Pando", "SCENICplus", "endpoint_summary.tsv", "endpoint_availability.tsv"))
  collisions <- paths[file.exists(paths)]
  if (length(collisions)) {
    stop("Archive existing formal results before a full rerun: ", paste(collisions, collapse = ", "))
  }
}
if (dir.exists(destination)) {
  stop("Formal GSE274113 control object already exists: ", destination)
}

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(Signac)
  library(GenomeInfoDb)
  library(rhdf5)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

required_inventory <- file.path(
  inventory_root, c("validation.tsv", "replicate_cell_alignment.tsv", "release_file_manifest.tsv")
)
fragment_exception_validation <- file.path(
  repo_root, "results", "gse274113", "inventory",
  "GSE274113_fragment_release_exception", "validation.tsv"
)
required_inventory <- c(required_inventory, fragment_exception_validation)
rna_guide_h5 <- file.path(
  data_root, paste0("GSE274113_filtered_feature_bc_matrix_", replicates, ".h5")
)
multiome_h5 <- file.path(
  data_root, paste0("GSE274113_rep", replicates, "_filtered_feature_bc_matrix.h5")
)
metadata_path <- file.path(data_root, "GSE274113_annotated_metadata.csv.gz")
control_timepoints <- c("day 7", "day 9", "day 11", "day 14")
required <- c(required_inventory, rna_guide_h5, multiome_h5, metadata_path)
if (any(!file.exists(required))) {
  stop("Complete GSE274113 inventory and RNA-ATAC matrices are required.")
}

inventory_validation <- data.table::fread(required_inventory[[1]])
if (!identical(colnames(inventory_validation), c("field", "value")) ||
  anyDuplicated(inventory_validation$field)) {
  stop("GSE274113 inventory validation schema is invalid.")
}
inventory_value <- function(field) {
  value <- inventory_validation$value[inventory_validation$field == field]
  if (length(value) != 1L) stop("Missing inventory field: ", field)
  value[[1]]
}
expected_peak_rows <- as.numeric(inventory_value("released_h5_peak_rows"))
expected_unique_peaks <- as.numeric(inventory_value("unique_released_h5_peak_coordinates"))
expected_consensus_all <- as.numeric(
  inventory_value("overlap_or_adjacent_merged_regions_all_seqlevels")
)
expected_consensus_standard <- as.numeric(
  inventory_value("overlap_or_adjacent_merged_regions_ucsc_standard")
)

find_one <- function(pattern) {
  paths <- list.files(
    extracted_root,
    pattern = pattern, full.names = TRUE, recursive = TRUE
  )
  if (length(paths) != 1L) stop("Expected one extracted file matching ", pattern)
  paths[[1]]
}

peak_annotation_paths <- vapply(
  replicates,
  function(rep) find_one(paste0("rep", rep, "_atac_peak_annotation\\.tsv\\.gz$")),
  character(1)
)

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
metadata <- metadata[, -cell_column_index, drop = FALSE]
if (anyNA(colnames(metadata)) || any(!nzchar(colnames(metadata))) ||
  anyDuplicated(colnames(metadata))) {
  stop("Released metadata has missing, empty or duplicated column names after cell-ID normalization.")
}
if (nrow(metadata) != 137604L || anyDuplicated(metadata$cell_id) ||
  !setequal(unique(metadata$replicate), paste0("rep", replicates)) ||
  !setequal(unique(metadata$Timepoint), control_timepoints)) {
  stop("Released GSE274113 metadata failed its full-data contract.")
}
control_metadata <- metadata[metadata$target == "NT", , drop = FALSE]
heldout_metadata <- metadata[metadata$target != "NT", , drop = FALSE]
if (nrow(control_metadata) != 11447L || nrow(heldout_metadata) != 126157L) {
  stop("Control/perturbation split differs from the complete release.")
}

h5_peak_names <- lapply(multiome_h5, function(path) {
  feature_names <- as.character(rhdf5::h5read(path, "matrix/features/name"))
  feature_types <- as.character(rhdf5::h5read(path, "matrix/features/feature_type"))
  peaks <- feature_names[feature_types == "Peaks"]
  if (!length(peaks) || anyDuplicated(peaks)) {
    stop("Released RNA+ATAC H5 has missing or duplicated peak names: ", path)
  }
  peaks
})
raw_peak_rows <- sum(vapply(h5_peak_names, length, integer(1)))
peak_keys <- unlist(h5_peak_names, use.names = FALSE)
unique_peak_coordinates <- length(unique(peak_keys))
all_peak_ranges <- Signac::StringToGRanges(
  regions = unique(peak_keys), sep = c(":", "-")
)
standard_seqlevels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
local_hg38_seqinfo <- GenomeInfoDb::seqinfo(
  BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)[standard_seqlevels]
if (!identical(GenomeInfoDb::seqlevels(local_hg38_seqinfo), standard_seqlevels) ||
  anyNA(GenomeInfoDb::seqlengths(local_hg38_seqinfo)) ||
  any(GenomeInfoDb::seqlengths(local_hg38_seqinfo) <= 0) ||
  any(GenomeInfoDb::genome(local_hg38_seqinfo) != "hg38")) {
  stop("The installed local hg38 Seqinfo failed validation.")
}
# Cross-replicate consensus ranges (reduce + standard seqlevels + hg38 Seqinfo)
# live in the package; the release audit below stays dataset-specific.
consensus_all <- multiCSN::consensus_peak_set(all_peak_ranges)
consensus_peaks <- multiCSN::consensus_peak_set(
  all_peak_ranges, seqlevels = standard_seqlevels, seqinfo = local_hg38_seqinfo
)
consensus_peak_names <- names(consensus_peaks)
if (raw_peak_rows != expected_peak_rows ||
  unique_peak_coordinates != expected_unique_peaks ||
  length(consensus_all) != expected_consensus_all ||
  length(consensus_peaks) != expected_consensus_standard ||
  anyDuplicated(consensus_peak_names)) {
  stop("Cross-replicate consensus peaks differ from the complete release audit.")
}
peak_annotation_counts <- vapply(
  peak_annotation_paths,
  function(path) nrow(data.table::fread(path, select = "chrom")),
  integer(1)
)
rm(all_peak_ranges, peak_keys)
gc(verbose = FALSE)

annotation <- Signac::GetGRangesFromEnsDb(
  ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
)
primary_seqlevel_map <- stats::setNames(
  c(paste0("chr", 1:22), "chrX", "chrY", "chrM"),
  c(as.character(1:22), "X", "Y", "MT")
)
primary_seqlevel_map <- primary_seqlevel_map[
  names(primary_seqlevel_map) %in% GenomeInfoDb::seqlevels(annotation)
]
annotation <- GenomeInfoDb::renameSeqlevels(annotation, primary_seqlevel_map)
GenomeInfoDb::genome(annotation) <- "hg38"

dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
stage <- tempfile(pattern = ".staging_GSE274113_control_full_", tmpdir = output_parent)
dir.create(stage)
rna_matrices <- vector("list", length(replicates))
atac_matrices <- vector("list", length(replicates))
metadata_rows <- vector("list", length(replicates))
validation_rows <- vector("list", length(replicates))
reference_rna_features <- NULL

cat("started_utc\t", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
for (index in seq_along(replicates)) {
  rep <- replicates[[index]]
  rep_name <- paste0("rep", rep)
  cat("replicate_started\t", rep_name, "\t", format(Sys.time(), usetz = TRUE), "\n", sep = "")
  counts <- Seurat::Read10X_h5(rna_guide_h5[[index]], use.names = TRUE, unique.features = TRUE)
  if (!is.list(counts) || !all(c("Gene Expression", "CRISPR Guide Capture") %in% names(counts))) {
    stop("Missing RNA or guide feature type for ", rep_name)
  }
  rna <- counts[["Gene Expression"]]
  rna_guide_h5_barcodes <- ncol(rna)
  rep_metadata <- control_metadata[control_metadata$replicate == rep_name, , drop = FALSE]
  annotated_replicate_cells <- sum(metadata$replicate == rep_name)
  rep_metadata$raw_barcode <- unname(sub(
    paste0("^", rep_name, "_"), "", rep_metadata$cell_id
  ))
  if (anyDuplicated(rep_metadata$raw_barcode) ||
    !all(rep_metadata$raw_barcode %in% colnames(rna))) {
    stop("Control metadata does not map exactly to ", rep_name, " RNA matrix.")
  }
  rna <- rna[, rep_metadata$raw_barcode, drop = FALSE]
  if (is.null(reference_rna_features)) {
    reference_rna_features <- rownames(rna)
  } else if (!identical(rownames(rna), reference_rna_features)) {
    stop("Released RNA feature order differs across replicates at ", rep_name)
  }
  if (!all(as.numeric(Matrix::colSums(rna)) == as.numeric(rep_metadata$nCount_RNA))) {
    stop("Released RNA counts disagree with metadata for ", rep_name)
  }
  multiome_counts <- Seurat::Read10X_h5(
    multiome_h5[[index]],
    use.names = TRUE, unique.features = TRUE
  )
  if (!is.list(multiome_counts) ||
    !all(c("Gene Expression", "Peaks") %in% names(multiome_counts))) {
    stop("Released RNA+ATAC H5 lacks a peak matrix for ", rep_name)
  }
  release_atac <- multiome_counts[["Peaks"]]
  rna_atac_h5_barcodes <- ncol(release_atac)
  if (!all(rep_metadata$raw_barcode %in% colnames(release_atac))) {
    stop("At least one control metadata cell is absent from ", rep_name, " ATAC matrix.")
  }
  release_atac <- release_atac[, rep_metadata$raw_barcode, drop = FALSE]
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
  nonstandard_peak_rows <- sum(!standard_peaks)
  nonstandard_counts <- sum(release_atac[!standard_peaks, , drop = FALSE])
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
  atac <- multiCSN::project_counts_to_consensus(
    release_atac, rownames(release_atac), consensus_peaks, type = "within"
  )
  if (!identical(rownames(atac), consensus_peak_names) ||
    !identical(colnames(atac), rep_metadata$raw_barcode)) {
    stop("Consensus H5 peak aggregation changed features or cells for ", rep_name)
  }
  colnames(rna) <- rep_metadata$cell_id
  colnames(atac) <- rep_metadata$cell_id
  rownames(rep_metadata) <- rep_metadata$cell_id
  validation_rows[[index]] <- data.frame(
    replicate = rep_name,
    timepoint = unique(rep_metadata$Timepoint),
    annotated_replicate_cells = annotated_replicate_cells,
    control_cells = ncol(rna),
    rna_atac_h5_barcodes = rna_atac_h5_barcodes,
    rna_atac_h5_barcodes_outside_annotation =
      rna_atac_h5_barcodes - annotated_replicate_cells,
    rna_guide_h5_barcodes = rna_guide_h5_barcodes,
    rna_guide_h5_barcodes_outside_annotation =
      rna_guide_h5_barcodes - annotated_replicate_cells,
    rna_features = nrow(rna),
    atac_features = nrow(atac),
    released_h5_replicate_peaks = length(h5_peak_names[[index]]),
    native_nonstandard_peaks_excluded = nonstandard_peak_rows,
    native_nonstandard_counts_excluded = nonstandard_counts,
    released_peak_annotation_rows = peak_annotation_counts[[index]],
    rna_nonzero = Matrix::nnzero(rna),
    atac_nonzero = Matrix::nnzero(atac),
    fragments_attached = FALSE,
    atac_quantification_source = "released_rna_atac_h5_consensus_aggregation",
    cell_cap = "none",
    feature_cap = "none",
    stringsAsFactors = FALSE
  )
  rna_matrices[[index]] <- rna
  atac_matrices[[index]] <- atac
  metadata_rows[[index]] <- rep_metadata
  rm(
    rna, atac, counts, multiome_counts, release_atac, release_ranges,
    within_replicate_overlaps
  )
  gc(verbose = FALSE)
  cat("replicate_completed\t", rep_name, "\t", format(Sys.time(), usetz = TRUE), "\n", sep = "")
}

rna_counts <- do.call(cbind, rna_matrices)
atac_counts <- do.call(cbind, atac_matrices)
combined_metadata <- do.call(rbind, metadata_rows)
combined_metadata <- combined_metadata[colnames(rna_counts), , drop = FALSE]
if (!identical(colnames(rna_counts), colnames(atac_counts)) ||
  !identical(colnames(rna_counts), rownames(combined_metadata)) ||
  anyDuplicated(colnames(rna_counts))) {
  stop("Sparse count assembly changed control-cell identity or order.")
}
rm(rna_matrices, atac_matrices, metadata_rows)
gc(verbose = FALSE)
combined <- Seurat::CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  meta.data = combined_metadata,
  min.cells = 0,
  min.features = 0,
  project = "GSE274113_all_NT_controls"
)
combined <- Seurat::NormalizeData(combined, assay = "RNA", verbose = FALSE)
combined[["ATAC"]] <- Signac::CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = local_hg38_seqinfo,
  annotation = annotation,
  min.cells = 0,
  min.features = 0
)
combined <- Signac::RunTFIDF(combined, assay = "ATAC")
if (ncol(combined) != 11447L ||
  !setequal(unique(combined$Timepoint), control_timepoints) ||
  !setequal(unique(combined$replicate), paste0("rep", replicates))) {
  stop("Combined GSE274113 control object lost cells, timepoints or replicates.")
}
combined@misc$full_multiome_20260803 <- list(
  dataset = "GSE274113 perturb-Multiome controls",
  scope = "all 11,447 NT controls across all 14 released replicates and four timepoints",
  heldout_scope = "all 126,157 non-NT cells represented by author perturbation-effect tables",
  split_design = "control-only network inference; perturbations held out for causal validation",
  cell_cap = "none",
  feature_cap = "none",
  atac_quantification_source = paste(
    "all 14 released RNA+ATAC H5 peak-count matrices aggregated to",
    "overlap-merged consensus regions"
  ),
  fragment_quantitative_input = "none; official rep2 fragment is truncated",
  genome = "hg38"
)
stored_rna_counts <- SeuratObject::LayerData(combined, assay = "RNA", layer = "counts")
stored_atac_counts <- SeuratObject::LayerData(combined, assay = "ATAC", layer = "counts")
expected_atac_ranges <- Signac::StringToGRanges(
  rownames(atac_counts),
  sep = c(":", "-")
)
stored_atac_ranges <- Signac::StringToGRanges(
  rownames(stored_atac_counts),
  sep = c(":", "-")
)
atac_coordinates_identical <-
  identical(
    as.character(GenomeInfoDb::seqnames(expected_atac_ranges)),
    as.character(GenomeInfoDb::seqnames(stored_atac_ranges))
  ) &&
    identical(
      GenomicRanges::start(expected_atac_ranges),
      GenomicRanges::start(stored_atac_ranges)
    ) &&
    identical(
      GenomicRanges::end(expected_atac_ranges),
      GenomicRanges::end(stored_atac_ranges)
    )
expected_atac_counts <- atac_counts
rownames(expected_atac_counts) <- rownames(stored_atac_counts)
if (!identical(rna_counts, stored_rna_counts) ||
  !isTRUE(atac_coordinates_identical) ||
  !identical(expected_atac_counts, stored_atac_counts)) {
  stop("Seurat construction changed the assembled sparse count matrices or coordinates.")
}
atac_counts <- stored_atac_counts
rm(
  stored_rna_counts, stored_atac_counts, expected_atac_counts,
  expected_atac_ranges, stored_atac_ranges
)
gc(verbose = FALSE)
object_path <- file.path(stage, "GSE274113_all_NT_controls_full_seurat.rds")
saveRDS(combined, object_path, compress = FALSE)
gene_ranges <- ensembldb::genes(
  EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
  return.type = "GRanges"
)
gene_seqlevel_map <- primary_seqlevel_map[
  names(primary_seqlevel_map) %in% GenomeInfoDb::seqlevels(gene_ranges)
]
gene_ranges <- GenomeInfoDb::renameSeqlevels(gene_ranges, gene_seqlevel_map)
gene_seqlevels <- standard_seqlevels[
  standard_seqlevels %in% GenomeInfoDb::seqlevels(gene_ranges)
]
gene_ranges <- GenomeInfoDb::keepSeqlevels(
  gene_ranges, gene_seqlevels,
  pruning.mode = "coarse"
)
GenomeInfoDb::seqinfo(gene_ranges) <- local_hg38_seqinfo[gene_seqlevels]
promoter_ranges <- GenomicRanges::promoters(
  gene_ranges,
  upstream = 2000L, downstream = 500L
)
accessible_consensus_peaks <- consensus_peaks[Matrix::rowSums(atac_counts) > 0]
accessible_consensus_peak_names <- Signac::GRangesToString(
  grange = accessible_consensus_peaks,
  sep = c(":", "-")
)
write_tsv(
  data.frame(peak = accessible_consensus_peak_names, stringsAsFactors = FALSE),
  file.path(stage, "accessible_consensus_peaks.tsv"), writer = "table"
)
accessibility_by_timepoint <- lapply(control_timepoints, function(timepoint) {
  timepoint_cells <- colnames(combined)[combined$Timepoint == timepoint]
  timepoint_accessible <- Matrix::rowSums(
    atac_counts[, timepoint_cells, drop = FALSE]
  ) > 0
  timepoint_slug <- gsub(" ", "_", timepoint, fixed = TRUE)
  write_tsv(
    data.frame(
      peak = consensus_peak_names[timepoint_accessible],
      stringsAsFactors = FALSE
    ),
    file.path(stage, paste0("accessible_consensus_peaks_", timepoint_slug, ".tsv")), writer = "table"
  )
  data.frame(
    timepoint = timepoint,
    control_cells = length(timepoint_cells),
    accessible_consensus_peaks = sum(timepoint_accessible),
    stringsAsFactors = FALSE
  )
})
accessibility_by_timepoint <- do.call(rbind, accessibility_by_timepoint)
if (any(accessibility_by_timepoint$accessible_consensus_peaks < 10000L)) {
  stop("A control timepoint has unexpectedly few accessible consensus peaks.")
}
write_tsv(
  accessibility_by_timepoint,
  file.path(stage, "accessible_consensus_peaks_by_timepoint.tsv"), writer = "table"
)
promoter_hits <- GenomicRanges::findOverlaps(
  promoter_ranges, accessible_consensus_peaks,
  ignore.strand = TRUE
)
supported_promoter_genes <- unique(as.character(
  S4Vectors::mcols(gene_ranges)$gene_name[S4Vectors::queryHits(promoter_hits)]
))
supported_promoter_genes <- sort(intersect(
  supported_promoter_genes[!is.na(supported_promoter_genes)],
  rownames(combined[["RNA"]])
))
if (length(supported_promoter_genes) < 10000L) {
  stop("Accessibility-supported promoter gene set is unexpectedly small.")
}
write_tsv(
  data.frame(gene = supported_promoter_genes, stringsAsFactors = FALSE),
  file.path(stage, "accessibility_supported_promoter_genes.tsv"), writer = "table"
)

write_sparse_group <- function(file, group, matrix, features, barcodes) {
  matrix <- methods::as(matrix, "dgCMatrix")
  rhdf5::h5createGroup(file, group)
  create_numeric <- function(name, values) {
    path <- paste0(group, "/", name)
    rhdf5::h5createDataset(
      file, path,
      dims = length(values),
      storage.mode = "integer",
      chunk = min(length(values), 1000000L),
      level = 6
    )
    rhdf5::h5write(as.integer(values), file, path)
  }
  create_numeric("data", matrix@x)
  create_numeric("indices", matrix@i)
  create_numeric("indptr", matrix@p)
  rhdf5::h5write(as.integer(dim(matrix)), file, paste0(group, "/shape"))
  rhdf5::h5write(as.character(features), file, paste0(group, "/features"))
  rhdf5::h5write(as.character(barcodes), file, paste0(group, "/barcodes"))
}
counts_h5 <- file.path(stage, "GSE274113_all_NT_controls_counts.h5")
rhdf5::h5createFile(counts_h5)
write_sparse_group(
  counts_h5, "rna",
  rna_counts,
  rownames(combined[["RNA"]]), colnames(combined)
)
write_sparse_group(
  counts_h5, "atac",
  atac_counts,
  rownames(combined[["ATAC"]]), colnames(combined)
)
rhdf5::H5close()
metadata_export <- combined[[]]
metadata_export$cell_id <- rownames(metadata_export)
data.table::fwrite(
  metadata_export, file.path(stage, "GSE274113_all_NT_controls_metadata.tsv.gz"),
  sep = "\t", quote = FALSE
)
validation <- do.call(rbind, validation_rows)
write_tsv(validation, file.path(stage, "replicate_object_validation.tsv"), writer = "table")
summary <- data.frame(
  field = c(
    "author_annotated_release_cells", "network_training_control_cells", "heldout_perturbation_cells",
    "rna_atac_h5_barcodes_total", "rna_atac_h5_barcodes_outside_annotation",
    "rna_guide_h5_barcodes_total", "rna_guide_h5_barcodes_outside_annotation",
    "replicates", "timepoints", "rna_release_features", "atac_consensus_features",
    "released_h5_peak_rows", "unique_released_h5_peak_coordinates",
    "released_peak_annotation_rows",
    "merged_peak_regions_all_seqlevels", "consensus_sequence_policy",
    "nonstandard_native_peak_rows_excluded", "nonstandard_native_counts_excluded",
    "accessible_consensus_peaks",
    "accessibility_supported_promoter_genes",
    "fragment_objects", "atac_quantification_source", "fragment_quantitative_input",
    "training_policy", "heldout_policy", "cell_cap",
    "feature_cap", "object_bytes", "counts_h5_bytes"
  ),
  value = c(
    nrow(metadata), ncol(combined), nrow(heldout_metadata),
    inventory_value("rna_atac_h5_barcodes_total"),
    inventory_value("multiome_barcodes_outside_annotated_metadata"),
    inventory_value("rna_guide_h5_barcodes_total"),
    inventory_value("guide_barcodes_outside_annotated_metadata"),
    length(replicates),
    length(unique(combined$Timepoint)), nrow(combined[["RNA"]]), nrow(combined[["ATAC"]]),
    raw_peak_rows, unique_peak_coordinates, sum(peak_annotation_counts),
    length(consensus_all),
    "overlap-merged release peaks on UCSC standard chromosomes; no numeric cap",
    sum(validation$native_nonstandard_peaks_excluded),
    sum(validation$native_nonstandard_counts_excluded),
    length(accessible_consensus_peaks),
    length(supported_promoter_genes),
    0L,
    "all 14 released RNA+ATAC H5 peak-count matrices aggregated to consensus",
    "none; official rep2 fragment is truncated", "all NT controls",
    "all non-NT perturbations", "none", "none", file.info(object_path)$size,
    file.info(counts_h5)$size
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary, file.path(stage, "control_object_validation.tsv"), writer = "table")
write_tsv(
  data.frame(
    file = c(
      metadata_path, rna_guide_h5, multiome_h5, peak_annotation_paths,
      required_inventory
    ),
    sha256 = unname(sha256_files(c(
      metadata_path, rna_guide_h5, multiome_h5, peak_annotation_paths,
      required_inventory
    ))),
    stringsAsFactors = FALSE
  ),
  file.path(stage, "input_SHA256SUMS.tsv"), writer = "table"
)
if (!file.rename(stage, destination)) stop("Unable to promote GSE274113 control object.")
cat("destination\t", destination, "\n", sep = "")
cat("completed_utc\t", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
print(summary)
