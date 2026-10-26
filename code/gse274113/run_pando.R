#!/usr/bin/env Rscript
run_pando <- function(cores = 1L, output_root = NULL) {

source("code/functions/network_io.R")

source("code/functions/runtime.R")
configure_task_runtime(comparator_library = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(Signac)
  library(Pando)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
data.table::setDTthreads(1L)

# ---- Parameters ----
options <- list(cores = cores, output_root = output_root)
cores <- options$cores
# ---------------------------------------------------------------
if (cores > 1L) {
  suppressPackageStartupMessages(library(foreach))
  doParallel::registerDoParallel(cores = cores)
  on.exit({doParallel::stopImplicitCluster(); foreach::registerDoSEQ()}, add = TRUE)
}
repo_root <- normalizePath(getwd(), mustWork = TRUE)
input_dir <- file.path(
  repo_root, "data", "gse274113", "objects", "GSE274113_control_full"
)
input_path <- file.path(input_dir, "GSE274113_all_NT_controls_full_seurat.rds")
input_validation <- file.path(input_dir, "control_object_validation.tsv")
destination <- if (is.null(options$output_root)) {
  file.path(repo_root, "runs", "gse274113", "method_fits", "Pando")
} else {
  normalizePath(options$output_root, mustWork = FALSE)
}
output_parent <- dirname(destination)
if (!all(file.exists(c(input_path, input_validation)))) {
  stop("GSE274113 control object and candidate data are required.")
}
if (dir.exists(destination)) stop("Formal Pando result already exists: ", destination)

timepoints <- c("day 7", "day 9", "day 11", "day 14")
pando_p_thresh <- 0.05
pando_rsq_thresh <- 0.1
pando_nvar_thresh <- 10L
pando_min_genes_per_module <- 5L
dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
stage <- tempfile(pattern = ".staging_Pando_GSE274113_", tmpdir = output_parent)
dir.create(stage)

object <- readRDS(input_path)
observed_timepoint_cells <- table(object$Timepoint)
if (!ncol(object) || anyDuplicated(colnames(object)) || anyNA(object$Timepoint) ||
  !setequal(names(observed_timepoint_cells), timepoints) ||
  !"target" %in% names(object[[]]) || anyNA(object$target) || any(object$target != "NT")) {
  stop("GSE274113 control object failed Pando input validation.")
}

data("phastConsElements20Mammals.UCSC.hg38", package = "Pando", envir = environment())
data("motifs", package = "Pando", envir = environment())
data("motif2tf", package = "Pando", envir = environment())

DefaultAssay(object) <- "RNA"
template <- Pando::initiate_grn(
  object,
  regions = phastConsElements20Mammals.UCSC.hg38,
  peak_assay = "ATAC",
  rna_assay = "RNA",
  exclude_exons = TRUE
)
template <- Pando::find_motifs(
  template,
  pfm = motifs,
  motif_tfs = motif2tf,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  verbose = TRUE
)
summary_rows <- vector("list", length(timepoints))
names(summary_rows) <- timepoints

for (index in seq_along(timepoints)) {
  timepoint <- timepoints[[index]]
  slug <- gsub(" ", "_", timepoint, fixed = TRUE)
  network_name <- paste0(slug, "_glm_network")
  time_dir <- file.path(stage, slug)
  dir.create(time_dir)
  cells <- colnames(template@data)[template@data$Timepoint == timepoint]
  if (!setequal(cells, colnames(object)[object$Timepoint == timepoint])) {
    stop("Unexpected control-cell count for ", timepoint)
  }
  fitted <- template
  fitted@data <- subset(template@data, cells = cells)
  DefaultAssay(fitted@data) <- "RNA"
  fitted@data <- Seurat::FindVariableFeatures(
    fitted@data,
    assay = "RNA",
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  selected_genes <- Seurat::VariableFeatures(fitted@data, assay = "RNA")
  if (!length(selected_genes)) stop("Pando variable-feature selection returned no genes.")
  fit_started <- Sys.time()
  fitted <- Pando::infer_grn(
    fitted,
    genes = selected_genes,
    network_name = network_name,
    peak_to_gene_method = "GREAT",
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    parallel = cores > 1L,
    tf_cor = 0.1,
    peak_cor = 0,
    method = "glm",
    family = "gaussian",
    interaction_term = ":",
    adjust_method = "fdr",
    scale = FALSE,
    verbose = TRUE
  )
  fit_seconds <- as.numeric(difftime(Sys.time(), fit_started, units = "secs"))
  coefficients <- as.data.frame(stats::coef(fitted))
  goodness_of_fit <- as.data.frame(Pando::gof(fitted))
  fitted <- Pando::find_modules(
    fitted,
    network = network_name,
    p_thresh = pando_p_thresh,
    rsq_thresh = pando_rsq_thresh,
    nvar_thresh = pando_nvar_thresh,
    min_genes_per_module = pando_min_genes_per_module
  )
  module_meta <- as.data.frame(Pando::NetworkModules(fitted)@meta)
  required_coefficient_columns <- c(
    "tf", "region", "target", "estimate", "pval", "padj"
  )
  required_fit_columns <- c("target", "rsq", "nvariables")
  required_module_columns <- c("tf", "target", "estimate", "pval", "padj")
  if (!all(required_coefficient_columns %in% names(coefficients)) ||
    !all(required_fit_columns %in% names(goodness_of_fit)) ||
    !all(required_module_columns %in% names(module_meta))) {
    stop("Pando coefficient, fit, or module schema changed at ", timepoint)
  }
  eligible_targets <- unique(as.character(goodness_of_fit$target[
    is.finite(goodness_of_fit$rsq) & goodness_of_fit$rsq > pando_rsq_thresh &
      is.finite(goodness_of_fit$nvariables) &
      goodness_of_fit$nvariables > pando_nvar_thresh
  ]))
  significant_coefficients <- data.table::as.data.table(coefficients)[
    target %chin% eligible_targets & (is.na(padj) | padj < pando_p_thresh)
  ]
  if (anyDuplicated(significant_coefficients[, .(tf, region, target)])) {
    stop(
      "Pando retained interaction terms have duplicated triplet keys at ",
      timepoint
    )
  }
  reconstructed_modules <- significant_coefficients[, .(
    estimate = sum(estimate), pval = min(pval), padj = min(padj)
  ), by = .(tf, target)]
  observed_modules <- data.table::as.data.table(module_meta)[
    , .(
      tf = as.character(tf), target = as.character(target),
      estimate = as.numeric(estimate), pval = as.numeric(pval),
      padj = as.numeric(padj)
    )
  ]
  data.table::setorder(reconstructed_modules, tf, target)
  data.table::setorder(observed_modules, tf, target)
  module_reconstruction_pass <- identical(
    reconstructed_modules[, .(tf, target)], observed_modules[, .(tf, target)]
  ) && isTRUE(all.equal(
    reconstructed_modules[, .(estimate, pval, padj)],
    observed_modules[, .(estimate, pval, padj)],
    tolerance = 1e-12, check.attributes = FALSE
  ))
  if (!module_reconstruction_pass) {
    stop(
      "Pando module metadata could not be reconstructed from coefficients and fit gates at ",
      timepoint
    )
  }
  write_tsv(data.frame(gene = selected_genes), file.path(time_dir, "native_variable_features.tsv"), na = "")
  write_tsv(coefficients, file.path(time_dir, "native_coefficients.tsv.gz"), na = "")
  write_tsv(
    significant_coefficients[, .(tf, region, target, estimate, pval, padj)],
    file.path(time_dir, "retained_native_triplets.tsv.gz"), na = ""
  )
  write_tsv(
    unique(significant_coefficients[, .(tf, region)])[order(tf, region)],
    file.path(time_dir, "retained_tf_region_binary.tsv.gz"), na = ""
  )
  write_tsv(
    unique(significant_coefficients[, .(region, target)])[order(region, target)],
    file.path(time_dir, "retained_region_gene_binary.tsv.gz"), na = ""
  )
  write_tsv(goodness_of_fit, file.path(time_dir, "native_goodness_of_fit.tsv"), na = "")
  write_tsv(module_meta, file.path(time_dir, "native_module_metadata.tsv"), na = "")
  saveRDS(Pando::GetNetwork(fitted, network = network_name), file.path(time_dir, "native_network_object.rds"), compress = FALSE)
  summary_rows[[index]] <- data.frame(
    timepoint = timepoint,
    control_cells = length(cells),
    release_rna_features = nrow(fitted@data[["RNA"]]),
    release_atac_features = nrow(fitted@data[["ATAC"]]),
    native_variable_features = length(selected_genes),
    genes_with_fit = if ("target" %in% names(goodness_of_fit)) length(unique(goodness_of_fit$target)) else nrow(goodness_of_fit),
    native_coefficients = nrow(coefficients),
    native_module_edges = nrow(module_meta),
    fit_seconds = fit_seconds,
    score_scope = "within_timepoint_native_endpoint_only",
    stringsAsFactors = FALSE
  )
  rm(
    fitted, coefficients, goodness_of_fit, module_meta,
    significant_coefficients, reconstructed_modules, observed_modules
  )
  gc(verbose = FALSE)
}
summary <- do.call(rbind, summary_rows)
if (!identical(summary$timepoint, timepoints) || sum(summary$control_cells) != ncol(object)) {
  stop("Pando timepoint suite failed control-cell accounting.")
}
write_tsv(summary, file.path(stage, "timepoint_run_summary.tsv"), na = "")
if (!file.rename(stage, destination)) stop("Unable to promote formal Pando result.")
print(summary)
}



if (sys.nframe() == 0L) run_pando()
