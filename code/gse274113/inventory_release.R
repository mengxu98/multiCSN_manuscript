source("code/functions/network_io.R")

inventory_release <- function(data_root = "data/gse274113/GSE274113",
    destination = "results/gse274113/inventory") {
  suppressPackageStartupMessages({library(data.table); library(rhdf5); library(Signac); library(multiCSN)})
  if (dir.exists(destination)) stop("Release inventory already exists: ", destination)
  metadata_path <- file.path(data_root, "GSE274113_annotated_metadata.csv.gz")
  metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
  column <- which(is.na(names(metadata)) | !nzchar(names(metadata)) | grepl("^Unnamed:", names(metadata)))
  if (length(column) != 1L || !all(c("replicate", "target", "Timepoint") %in% names(metadata))) stop("Invalid released metadata schema.")
  ids <- as.character(metadata[[column]])
  metadata <- metadata[!is.na(ids) & nzchar(trimws(ids)), , drop = FALSE]
  metadata$cell_id <- as.character(metadata[[column]])
  if (!nrow(metadata) || anyDuplicated(metadata$cell_id) ||
      anyNA(metadata$target) || any(!nzchar(metadata$target)) ||
      !any(metadata$target == "NT") || !any(metadata$target != "NT") ||
      anyNA(metadata$replicate) ||
      any(!grepl("^rep[1-9][0-9]*$", metadata$replicate))) stop("Released cell universe or control split changed.")
  replicates <- sort(unique(as.integer(sub("^rep", "", metadata$replicate))))
  rna_paths <- file.path(data_root, paste0("GSE274113_filtered_feature_bc_matrix_", replicates, ".h5"))
  atac_paths <- file.path(data_root, paste0("GSE274113_rep", replicates, "_filtered_feature_bc_matrix.h5"))
  extracted <- list.files(file.path(data_root, "extracted"), recursive = TRUE, full.names = TRUE)
  annotation_paths <- vapply(replicates, function(replicate) {
    hits <- extracted[grepl(paste0("rep", replicate, "_atac_peak_annotation\\.tsv\\.gz$"), extracted)]
    if (length(hits) != 1L) stop("Expected one peak annotation for replicate ", replicate)
    hits
  }, character(1))
  supplements <- file.path(data_root, "supplementary", c("NIHMS2076308-supplement-Table_S2.xlsx", "NIHMS2076308-supplement-Table_S3.xlsx"))
  inputs <- c(metadata_path, rna_paths, atac_paths, annotation_paths, supplements)
  if (any(!file.exists(inputs))) stop("Incomplete raw GSE274113 release.")
  alignment <- vector("list", length(replicates)); peaks <- vector("list", length(replicates))
  for (index in seq_along(replicates)) {
    replicate <- paste0("rep", replicates[[index]])
    rna <- as.character(h5read(rna_paths[[index]], "matrix/barcodes"))
    atac <- as.character(h5read(atac_paths[[index]], "matrix/barcodes"))
    features <- as.character(h5read(atac_paths[[index]], "matrix/features/name"))
    types <- as.character(h5read(atac_paths[[index]], "matrix/features/feature_type"))
    if (anyDuplicated(rna) || anyDuplicated(atac) || length(features) != length(types)) stop("Invalid released H5 barcode/feature schema.")
    peaks[[index]] <- features[types == "Peaks"]
    if (!length(peaks[[index]]) || anyDuplicated(peaks[[index]])) stop("Invalid released ATAC peak universe.")
    current <- metadata[metadata$replicate == replicate, , drop = FALSE]
    barcodes <- sub(paste0("^", replicate, "_"), "", current$cell_id)
    if (anyDuplicated(barcodes) || !all(barcodes %in% rna) || !all(barcodes %in% atac) || !all(rna %in% atac)) stop("Released metadata/RNA/ATAC cells do not align for ", replicate)
    annotation <- fread(annotation_paths[[index]], select = "chrom")
    if (!nrow(annotation) || anyNA(annotation$chrom)) stop("Invalid released peak annotation.")
    alignment[[index]] <- data.table(replicate = replicate, annotated_cells = nrow(current),
      control_cells = sum(current$target == "NT"), rna_guide_h5_barcodes = length(rna),
      rna_atac_h5_barcodes = length(atac), guide_barcodes_outside_annotated_metadata = sum(!rna %in% barcodes),
      multiome_barcodes_outside_annotated_metadata = sum(!atac %in% barcodes),
      released_h5_peak_rows = length(peaks[[index]]), peak_annotation_rows = nrow(annotation))
  }
  alignment <- rbindlist(alignment)
  coordinates <- unique(unlist(peaks, use.names = FALSE))
  ranges <- Signac::StringToGRanges(coordinates, sep = c(":", "-"))
  standard <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  seqinfo <- GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)[standard]
  values <- c(released_h5_peak_rows = sum(alignment$released_h5_peak_rows),
    unique_released_h5_peak_coordinates = length(coordinates),
    overlap_or_adjacent_merged_regions_all_seqlevels = length(multiCSN::consensus_peak_set(ranges)),
    overlap_or_adjacent_merged_regions_ucsc_standard = length(multiCSN::consensus_peak_set(ranges, seqlevels = standard, seqinfo = seqinfo)),
    rna_atac_h5_barcodes_total = sum(alignment$rna_atac_h5_barcodes),
    multiome_barcodes_outside_annotated_metadata = sum(alignment$multiome_barcodes_outside_annotated_metadata),
    rna_guide_h5_barcodes_total = sum(alignment$rna_guide_h5_barcodes),
    guide_barcodes_outside_annotated_metadata = sum(alignment$guide_barcodes_outside_annotated_metadata),
    annotated_metadata_cells = nrow(metadata), control_cells = sum(metadata$target == "NT"))
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(".release-inventory-", tmpdir = dirname(destination))
  dir.create(stage, recursive = TRUE)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  write_tsv(data.table(field = names(values), value = as.character(values)), file.path(stage, "validation.tsv"))
  write_tsv(alignment, file.path(stage, "replicate_cell_alignment.tsv"))
  if (!file.rename(stage, destination)) stop("Unable to publish the release inventory.")
  invisible(destination)
}

if (sys.nframe() == 0L) inventory_release()
