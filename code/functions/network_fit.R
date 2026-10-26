#!/usr/bin/env Rscript

source("code/functions/runtime.R")
configure_task_runtime()

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(Signac)
  library(multiCSN)
})
data.table::setDTthreads(1L)

source("code/functions/network_io.R")

# Shared fitting and publication contract for state and timepoint networks.
fit_state_networks <- function(
  dataset,
  input,
  output,
  state,
  expected_cells = NULL,
  response_chunk = 1024L,
  profile = "gse194122",
  cell_manifest = NULL,
  cores = 1L
) {
  arguments <- list(
    dataset = dataset, input = input, output = output,
    state = state, expected_cells = expected_cells, response_chunk = response_chunk,
    cell_manifest = cell_manifest, profile = profile, cores = cores
  )
  input_path <- normalizePath(arguments$input, mustWork = TRUE)
  output_root <- normalizePath(arguments$output, mustWork = FALSE)
  if (file.exists(output_root)) stop("Output exists and will not be overwritten: ", output_root)
  dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
  stage_root <- tempfile(pattern = paste0(".", basename(output_root), "-"), tmpdir = dirname(output_root))
  dir.create(stage_root)
  on.exit(unlink(stage_root, recursive = TRUE), add = TRUE)

  object <- readRDS(input_path)
  params <- multiCSN::Params(object)
  metadata <- object[[]]
  is_timepoint <- arguments$profile == "gse274113"
  if (is_timepoint && arguments$dataset == "GSE274113" &&
    (!identical(params$filter_mode, "unfiltered") ||
      !"target" %in% names(metadata) || anyNA(object$target) || any(object$target != "NT") || !setequal(
      unique(as.character(object$Timepoint)),
      c("day 7", "day 9", "day 11", "day 14")
    ))) {
    stop("Expected unfiltered NT controls across the declared timepoints.")
  }
  catalogue_type <- if (is_timepoint) arguments$state else "aggregate"
  state_column <- if (is_timepoint) "Timepoint" else "cell_type"
  catalogue_cells <- multiCSN::get_attribute(object, celltypes = catalogue_type, attribute = "cells")
  if (!state_column %in% colnames(metadata)) stop("State metadata is absent.")
  metadata_cells <- rownames(metadata)[which(as.character(metadata[[state_column]]) == arguments$state)]
  if (anyNA(metadata_cells) || anyDuplicated(catalogue_cells) ||
      !all(metadata_cells %in% catalogue_cells) ||
      (is_timepoint && !setequal(catalogue_cells, metadata_cells))) {
    stop("Candidate catalogue does not match the state metadata.")
  }
  cells <- intersect(catalogue_cells, metadata_cells)
  cell_manifest_path <- NULL
  if (!is.null(arguments$cell_manifest) && nzchar(arguments$cell_manifest)) {
    cell_manifest_path <- normalizePath(arguments$cell_manifest, mustWork = TRUE)
    requested_cells <- read_fit_checkpoint(cell_manifest_path)
    if (!identical(names(requested_cells), c("cell_id", "state")) ||
      is.null(arguments$expected_cells) ||
      nrow(requested_cells) != arguments$expected_cells ||
      anyDuplicated(requested_cells$cell_id) ||
      any(requested_cells$state != arguments$state)) {
      stop("Invalid cell manifest: ", arguments$state, ".")
    }
    if (any(!requested_cells$cell_id %in% cells)) {
      stop("Manifest cells are outside the candidate set.")
    }
    cells <- requested_cells$cell_id
  }
  targets_requested <- if (is_timepoint) {
    rownames(object[[params$rna_assay]])
  } else {
    multiCSN::get_attribute(object, celltypes = catalogue_type, attribute = "genes")
  }
  peaks_requested <- multiCSN::get_attribute(object, celltypes = catalogue_type, attribute = "peaks")
  regulators_requested <- multiCSN::get_attribute(object, attribute = "tfs")
  if (!length(cells) ||
    (!is.null(arguments$expected_cells) && length(cells) != arguments$expected_cells) || anyDuplicated(cells) ||
    !length(targets_requested) || !length(peaks_requested) || !length(regulators_requested)) {
    stop("Invalid candidate catalogue.")
  }
  # The layered fit itself (timepoint re-normalization, motif- and
  # domain-gated candidate construction, the three Greedy-l0 layers and the
  # factorized chain projection) is implemented by multiCSN::fit_layered_network().
  fit <- multiCSN::fit_layered_network(
    object,
    cell_group = catalogue_type,
    cells = cells,
    min_detected = if (is_timepoint) 20L else 0L,
    min_peak_cells = 3L,
    upstream = 100000,
    downstream = 0,
    tf_region_scope = if (is_timepoint) "motif_peaks" else "domain_peaks",
    exclude_response_alias = is_timepoint,
    renormalize = is_timepoint,
    response_chunk = arguments$response_chunk,
    cores = arguments$cores,
    sort_regulators = !is_timepoint,
    checkpoint_dir = stage_root,
    verbose = TRUE
  )

  multiCSN::write_network_endpoints(fit, stage_root, overwrite = TRUE)
  # The component tables above are the complete result. Per-target chunk checkpoints
  # only serve resume during fitting and are not part of the promoted output.
  part_dirs <- list.dirs(stage_root, recursive = FALSE, full.names = TRUE)
  part_dirs <- part_dirs[grepl("_parts$", part_dirs)]
  if (length(part_dirs)) unlink(part_dirs, recursive = TRUE, force = TRUE)
  if (!file.rename(stage_root, output_root)) stop("Unable to promote complete network result.")
  invisible(data.table(cell_id = cells, state = arguments$state))
}
