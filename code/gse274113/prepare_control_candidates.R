#!/usr/bin/env Rscript

source("code/functions/network_io.R")

source("code/functions/runtime.R")
configure_task_runtime()

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
  library(Signac)
  library(multiCSN)
  library(motifmatchr)
})

repo_root <- normalizePath(getwd(), mustWork = TRUE)
object_root <- file.path(
  repo_root, "data", "gse274113", "objects", "GSE274113_control_full"
)
input_path <- file.path(object_root, "GSE274113_all_NT_controls_full_seurat.rds")
input_validation <- file.path(object_root, "control_object_validation.tsv")
motif_dir <- file.path(repo_root, "data", "resources", "motifs")
destination <- file.path(object_root, "candidates")
required <- c(
  input_path,
  input_validation,
  file.path(motif_dir, c("motifs.rda", "motif2tf.rda"))
)
if (any(!file.exists(required))) stop("Complete GSE274113 control object and motif resources are required.")
if (dir.exists(destination)) stop("Formal GSE274113 candidates already exist: ", destination)

# find_motifs inherits these scanner defaults; reject a changed configuration.
scanner <- methods::selectMethod("matchMotifs", c("PWMatrixList", "GenomicRanges"))
scanner_defaults <- formals(eval(body(scanner)[[2L]][[3L]], environment(scanner)))
if (!identical(eval(scanner_defaults$p.cutoff), 5e-5) ||
    !identical(eval(scanner_defaults$bg)[[1L]], "subject") ||
    !identical(eval(scanner_defaults$w), 7)) {
  stop("motifmatchr defaults differ from the declared motif scan settings.")
}

stage <- tempfile(pattern = ".staging_candidates_", tmpdir = object_root)
dir.create(stage)
cat("started_utc\t", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
object <- readRDS(input_path)
timepoints <- c("day 7", "day 9", "day 11", "day 14")
if (ncol(object) != 11447L ||
  !setequal(unique(as.character(object$Timepoint)), timepoints) ||
  !identical(
    object@misc$full_multiome_20260803$split_design,
    "control-only network inference; perturbations held out for causal validation"
  )) {
  stop("GSE274113 control object failed split validation.")
}
rna_features <- rownames(object[["RNA"]])
atac_features <- rownames(object[["ATAC"]])
object <- multiCSN::initiate_object(
  object,
  group.by = "Timepoint",
  targets = rna_features,
  filter_mode = "unfiltered",
  filter_by = "aggregate",
  rna_assay = "RNA",
  peak_assay = "ATAC",
  regions = NULL,
  regions_extend = 0,
  exclude_exons = TRUE,
  verbose = TRUE
)
# Load the declared resources explicitly; find_motifs takes matrices and a map.
load_resource <- function(path) {
  resource <- new.env(parent = emptyenv())
  keys <- load(path, envir = resource)
  if (length(keys) != 1L) stop("Expected one resource object in ", path)
  resource[[keys[[1]]]]
}
object <- multiCSN::find_motifs(
  object,
  pfm = load_resource(file.path(motif_dir, "motifs.rda")),
  motif_tfs = load_resource(file.path(motif_dir, "motif2tf.rda")),
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
  backend = "motifmatchr",
  verbose = TRUE
)
observed_timepoints <- multiCSN::get_attribute(object, attribute = "celltypes")
selected_peaks <- unique(unlist(lapply(
  observed_timepoints,
  function(current) multiCSN::get_attribute(object, celltypes = current, attribute = "peaks")
), use.names = FALSE))
selected_genes <- unique(unlist(lapply(
  observed_timepoints,
  function(current) multiCSN::get_attribute(object, celltypes = current, attribute = "genes")
), use.names = FALSE))
tfs <- multiCSN::get_attribute(object, attribute = "tfs")
motif_matrix <- multiCSN::NetworkRegions(object)@motifs@data
genome_mapping <- multiCSN::map_peak_ranges_to_genome(
  multiCSN::parse_peak_ranges(selected_peaks),
  BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)
if (!setequal(observed_timepoints, timepoints) ||
  !identical(sort(selected_genes), sort(rna_features)) ||
  !all(selected_peaks %in% atac_features) ||
  !identical(sort(rownames(motif_matrix)), sort(selected_peaks)) ||
  !length(tfs)) {
  stop("GSE274113 timepoint candidate construction failed.")
}

candidate_path <- file.path(stage, "GSE274113_NT_timepoint_candidates.rds")
saveRDS(object, candidate_path, compress = FALSE)
timepoint_counts <- table(factor(object$Timepoint, levels = timepoints))
validation <- data.frame(
  field = c(
    "dataset", "network_training_policy", "heldout_policy", "author_annotated_release_cells",
    "training_control_cells", "heldout_perturbation_cells", "timepoints",
    paste0("training_cells_", gsub(" ", "_", timepoints)),
    "input_rna_features", "input_atac_features", "initialized_target_genes",
    "initialized_peaks", "peaks_removed_by_declared_region_rule", "candidate_tfs",
    "motifs", "motif_peak_pairs", "motif_nonzero_hits",
    "motif_genome_supported_peaks", "motif_genome_unsupported_peaks",
    "motif_seqname_remapped_peaks", "filter_mode",
    "grouping", "exclude_exons", "motif_method", "motif_score_mode", "cell_cap",
    "target_cap", "peak_cap", "candidate_edge_cap"
  ),
  value = c(
    "GSE274113", "all NT controls", "all non-NT perturbations", 137604, 11447,
    126157, length(timepoints), as.integer(timepoint_counts), length(rna_features),
    length(atac_features), length(selected_genes), length(selected_peaks),
    length(setdiff(atac_features, selected_peaks)), length(tfs), ncol(motif_matrix),
    as.double(nrow(motif_matrix)) * as.double(ncol(motif_matrix)),
    Matrix::nnzero(motif_matrix), sum(genome_mapping$supported),
    sum(!genome_mapping$supported),
    sum(genome_mapping$original_seqnames != genome_mapping$mapped_seqnames),
    "unfiltered", "Timepoint", TRUE, "motifmatchr", "binary",
    "none", "none", "none", "none"
  ),
  stringsAsFactors = FALSE
)
write_tsv(validation, file.path(stage, "candidate_validation.tsv"), writer = "table")

if (!file.rename(stage, destination)) stop("Unable to promote GSE274113 candidates.")
cat("destination\t", destination, "\n", sep = "")
cat("completed_utc\t", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
print(validation)
