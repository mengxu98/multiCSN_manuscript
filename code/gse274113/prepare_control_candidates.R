#!/usr/bin/env Rscript

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
motif_dir <- file.path(repo_root, "data", "resources", "motifs")
destination <- file.path(object_root, "candidates")
if (dir.exists(destination)) stop("Formal GSE274113 candidates already exist: ", destination)

stage <- tempfile(pattern = ".staging_candidates_", tmpdir = object_root)
dir.create(stage)
object <- readRDS(input_path)
timepoints <- c("day 7", "day 9", "day 11", "day 14")
if (!ncol(object) || anyDuplicated(colnames(object)) ||
  !setequal(unique(as.character(object$Timepoint)), timepoints) ||
  !"target" %in% colnames(object[[]]) || anyNA(object$target) || any(object$target != "NT")) {
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
if (!setequal(observed_timepoints, timepoints) ||
  !identical(sort(selected_genes), sort(rna_features)) ||
  !all(selected_peaks %in% atac_features) ||
  !identical(sort(rownames(motif_matrix)), sort(selected_peaks)) ||
  !length(tfs)) {
  stop("GSE274113 timepoint candidate construction failed.")
}

candidate_path <- file.path(stage, "GSE274113_NT_timepoint_candidates.rds")
saveRDS(object, candidate_path, compress = FALSE)
if (!file.rename(stage, destination)) stop("Unable to promote GSE274113 candidates.")
