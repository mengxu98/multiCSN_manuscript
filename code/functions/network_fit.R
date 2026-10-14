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
  validation,
  output,
  state,
  expected_cells,
  response_chunk = 1024L,
  profile = "gse194122",
  cell_manifest = NULL,
  cores = 1L
) {
  arguments <- list(
    dataset = dataset, input = input, validation = validation, output = output,
    state = state, expected_cells = expected_cells, response_chunk = response_chunk,
    cell_manifest = cell_manifest, profile = profile, cores = cores
  )
  script_path <- normalizePath("code/functions/network_fit.R", mustWork = TRUE)
  multicsn_package_root <- normalizePath(system.file(package = "multiCSN"), mustWork = TRUE)
  csn_package_root <- normalizePath(system.file(package = "inferCSN"), mustWork = TRUE)
  if (!"fit_layered_network" %in% getNamespaceExports("multiCSN")) {
    stop(
      "The installed multiCSN does not export fit_layered_network(). ",
      "Install the required packages using code/requirements/install_csn_packages.R."
    )
  }

  input_path <- normalizePath(arguments$input, mustWork = TRUE)
  validation_path <- normalizePath(arguments$validation, mustWork = TRUE)
  output_root <- normalizePath(arguments$output, mustWork = FALSE)
  if (file.exists(output_root)) stop("Output exists and will not be overwritten: ", output_root)
  dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
  stage_root <- paste0(output_root, ".staging")
  if (!dir.exists(stage_root)) dir.create(stage_root, recursive = FALSE)

  object <- readRDS(input_path)
  params <- multiCSN::Params(object)
  is_timepoint <- arguments$profile == "gse274113"
  if (is_timepoint && arguments$dataset == "GSE274113" &&
    (ncol(object) != 11447L || !setequal(
      unique(as.character(object$Timepoint)),
      c("day 7", "day 9", "day 11", "day 14")
    ))) {
    stop("Expected all 11,447 NT controls across the four declared timepoints.")
  }
  catalogue_type <- if (is_timepoint) arguments$state else "aggregate"
  state_column <- if (is_timepoint) "Timepoint" else "cell_type"
  catalogue_cells <- multiCSN::get_attribute(object, celltypes = catalogue_type, attribute = "cells")
  if (!state_column %in% colnames(object@meta.data)) stop("State metadata is absent.")
  cells <- intersect(
    catalogue_cells,
    rownames(object@meta.data)[as.character(object@meta.data[[state_column]]) == arguments$state]
  )
  cell_manifest_path <- NULL
  if (!is.null(arguments$cell_manifest) && nzchar(arguments$cell_manifest)) {
    cell_manifest_path <- normalizePath(arguments$cell_manifest, mustWork = TRUE)
    requested_cells <- read_fit_checkpoint(cell_manifest_path)
    if (!identical(names(requested_cells), c("cell_id", "state")) ||
      nrow(requested_cells) != arguments$expected_cells ||
      anyDuplicated(requested_cells$cell_id) ||
      any(requested_cells$state != arguments$state)) {
      stop("Exact cell manifest is invalid for state ", arguments$state, ".")
    }
    if (any(!requested_cells$cell_id %in% cells)) {
      stop("Exact cell manifest contains cells outside the complete state catalogue.")
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
  if (length(cells) != arguments$expected_cells || anyDuplicated(cells) ||
    !length(targets_requested) || !length(peaks_requested) || !length(regulators_requested)) {
    stop("Complete candidate catalogue validation failed.")
  }
  # The layered fit itself (timepoint re-normalization, motif- and
  # domain-gated candidate construction, the three Greedy-l0 layers and the
  # factorized chain projection) is implemented by multiCSN::fit_layered_network().
  input_files <- unique(normalizePath(c(
    input_path, validation_path, script_path,
    "code/functions/network_io.R", "code/functions/table_io.R",
    "code/functions/network_fit.R",
    cell_manifest_path,
    file.path(multicsn_package_root, c("DESCRIPTION", "NAMESPACE")),
    list.files(file.path(multicsn_package_root, "R"), full.names = TRUE, recursive = TRUE),
    list.files(file.path(multicsn_package_root, "libs"), full.names = TRUE, recursive = TRUE),
    file.path(csn_package_root, c("DESCRIPTION", "NAMESPACE")),
    list.files(file.path(csn_package_root, "R"), full.names = TRUE, recursive = TRUE),
    list.files(file.path(csn_package_root, "libs"), full.names = TRUE, recursive = TRUE)
  ), mustWork = TRUE))
  current_manifest <- data.table(
    file = input_files, bytes = file.info(input_files)$size,
    sha256 = unname(sha256_files(input_files))
  )
  settings <- arguments[setdiff(names(arguments), "output")]
  settings$input <- input_path
  settings$validation <- validation_path
  settings$cell_manifest <- cell_manifest_path
  identity <- list(settings = settings, inputs = current_manifest, cells = cells)
  manifest_path <- file.path(stage_root, "input_identity.rds")
  if (!file.exists(manifest_path) && length(list.files(stage_root, "_parts$"))) {
    stop("Staging checkpoints lack input identity; use a new output directory.")
  }
  if (file.exists(manifest_path)) {
    previous_manifest <- readRDS(manifest_path)
    if (!isTRUE(all.equal(previous_manifest, identity,
      check.attributes = FALSE
    ))) {
      stop("Staging input, code or package identity changed; use a new output directory.")
    }
  }
  saveRDS(identity, manifest_path)

  # The layered multi-omic fit itself lives in the multiCSN package; this script
  # only validates the catalogue, pins the input identity, materialises the
  # publication contract and promotes the staged result.
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

  tf_gene <- fit$tf_gene
  tf_region <- fit$tf_region
  region_gene <- fit$region_gene
  chains <- fit$triplets
  mediated <- fit$mediated_tf_gene
  # The publication contract names the TF-region response column `target`; the
  # package uses the semantic name `region`.
  tf_region_accounting <- fit$accounting$tf_region
  names(tf_region_accounting)[names(tf_region_accounting) == "region"] <- "target"

  write_tsv(tf_gene, file.path(stage_root, "tf_gene.tsv"))
  write_tsv(tf_region, file.path(stage_root, "tf_region.tsv"))
  write_tsv(region_gene, file.path(stage_root, "region_gene.tsv"))
  write_tsv(chains, file.path(stage_root, "triplets.tsv"))
  write_tsv(mediated, file.path(stage_root, "mediated_tf_gene.tsv"))
  write_tsv(fit$accounting$tf_gene, file.path(stage_root, "tf_gene_accounting.tsv"))
  write_tsv(tf_region_accounting, file.path(stage_root, "tf_region_accounting.tsv"))
  write_tsv(fit$accounting$region_gene, file.path(stage_root, "region_gene_accounting.tsv"))
  write_tsv(fit$summary, file.path(stage_root, "endpoint_summary.tsv"))
  # The component tables above are the complete result. Per-target chunk checkpoints
  # only serve resume during fitting and are not part of the promoted output.
  part_dirs <- list.dirs(stage_root, recursive = FALSE, full.names = TRUE)
  part_dirs <- part_dirs[grepl("_parts$", part_dirs)]
  if (length(part_dirs)) unlink(part_dirs, recursive = TRUE, force = TRUE)
  unlink(manifest_path)
  outputs <- sort(list.files(stage_root, recursive = TRUE, full.names = TRUE), method = "radix")
  outputs <- outputs[!file.info(outputs)$isdir]
  output_relative <- substring(outputs, nchar(stage_root) + 2L)
  identity$outputs <- data.table(file = output_relative, sha256 = unname(sha256_files(outputs)))
  saveRDS(identity, file.path(stage_root, "run_provenance.rds"))
  source("code/functions/network_contract.R")
  validate_network_result(stage_root, settings)
  if (!file.rename(stage_root, output_root)) stop("Unable to promote complete network result.")
  cat("network_fit_complete\tPASS\t", output_root, "\n", sep = "")
  invisible(data.table(cell_id = cells, state = arguments$state))
}
