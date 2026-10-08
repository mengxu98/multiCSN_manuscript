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
if (!identical(cores, 1L)) stop("This benchmark requires cores = 1.")
repo_root <- normalizePath(getwd(), mustWork = TRUE)
script_path <- normalizePath("code/gse274113/run_pando.R", mustWork = TRUE)
execution_root <- normalizePath(file.path(repo_root, ".."), mustWork = TRUE)
input_dir <- file.path(
  repo_root, "data", "gse274113", "objects", "GSE274113_control_full"
)
input_path <- file.path(input_dir, "GSE274113_all_NT_controls_full_seurat.rds")
input_validation <- file.path(input_dir, "control_object_validation.tsv")
source_commit_path <- file.path(
  repo_root, "logs", "gse274113", "comparator_install",
  "Pando_source_commit.txt"
)
destination <- if (is.null(options$output_root)) {
  file.path(repo_root, "runs", "gse274113", "method_fits", "Pando")
} else {
  normalizePath(options$output_root, mustWork = FALSE)
}
output_parent <- dirname(destination)
if (!all(file.exists(c(input_path, input_validation, source_commit_path)))) {
  stop("Complete GSE274113 control object and Pando provenance are required.")
}
if (dir.exists(destination)) stop("Formal Pando result already exists: ", destination)

safe_slot <- function(object, expression, fallback = data.frame()) {
  tryCatch(force(expression), error = function(error) fallback)
}

timepoints <- c("day 7", "day 9", "day 11", "day 14")
expected_cells <- c("day 7" = 3304L, "day 9" = 4167L, "day 11" = 1903L, "day 14" = 2073L)
pando_p_thresh <- 0.05
pando_rsq_thresh <- 0.1
pando_nvar_thresh <- 10L
pando_min_genes_per_module <- 5L
dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
stage <- tempfile(pattern = ".staging_Pando_GSE274113_", tmpdir = output_parent)
dir.create(stage)

cat("started_utc\t", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
input_sha256 <- unname(sha256_files(input_path))
object <- readRDS(input_path)
observed_timepoint_cells <- table(object$Timepoint)
if (ncol(object) != 11447L ||
  !setequal(names(observed_timepoint_cells), timepoints) ||
  !identical(as.integer(observed_timepoint_cells[timepoints]), unname(expected_cells))) {
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
template_regions <- Pando::NetworkRegions(template)
template_validation <- data.frame(
  field = c(
    "input_control_cells", "input_rna_features", "input_atac_features",
    "candidate_ranges_after_phastcons_and_exon_filter", "candidate_peak_mappings",
    "motif_rows", "motif_columns", "motif_nonzero", "motif_tf_rows",
    "motif_tf_columns", "motif_tf_nonzero"
  ),
  value = c(
    ncol(object), nrow(object[["RNA"]]), nrow(object[["ATAC"]]),
    length(template_regions@ranges), length(template_regions@peaks),
    nrow(template_regions@motifs@data), ncol(template_regions@motifs@data),
    Matrix::nnzero(template_regions@motifs@data), nrow(template_regions@tfs),
    ncol(template_regions@tfs), Matrix::nnzero(template_regions@tfs)
  ),
  stringsAsFactors = FALSE
)
write_tsv(template_validation, file.path(stage, "shared_candidate_validation.tsv"), na = "")

summary_rows <- vector("list", length(timepoints))
names(summary_rows) <- timepoints
reconstruction_rows <- vector("list", length(timepoints))
names(reconstruction_rows) <- timepoints

for (index in seq_along(timepoints)) {
  timepoint <- timepoints[[index]]
  slug <- gsub(" ", "_", timepoint, fixed = TRUE)
  network_name <- paste0(slug, "_glm_network")
  time_dir <- file.path(stage, slug)
  dir.create(time_dir)
  cells <- colnames(template@data)[template@data$Timepoint == timepoint]
  if (length(cells) != expected_cells[[timepoint]]) {
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
    parallel = FALSE,
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
  module_meta <- safe_slot(
    fitted,
    as.data.frame(Pando::NetworkModules(fitted)@meta)
  )
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
  reconstruction_rows[[index]] <- data.frame(
    timepoint = timepoint,
    fitted_targets = length(unique(goodness_of_fit$target)),
    eligible_targets_rsq_and_nvariables = length(eligible_targets),
    coefficients_before_module_gates = nrow(coefficients),
    coefficients_after_native_gates = nrow(significant_coefficients),
    reconstructed_tf_target_edges = nrow(reconstructed_modules),
    observed_tf_target_edges = nrow(observed_modules),
    reconstruction_tolerance = 1e-12,
    status = "PASS",
    stringsAsFactors = FALSE
  )
  write_tsv(
    reconstruction_rows[[index]],
    file.path(time_dir, "native_module_reconstruction.tsv"), na = ""
  )
  write_tsv(data.frame(gene = selected_genes), file.path(time_dir, "native_variable_features.tsv"), na = "")
  saveRDS(coefficients, file.path(time_dir, "native_coefficients.rds"), compress = FALSE)
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
  saveRDS(goodness_of_fit, file.path(time_dir, "native_goodness_of_fit.rds"), compress = FALSE)
  write_tsv(goodness_of_fit, file.path(time_dir, "native_goodness_of_fit.tsv"), na = "")
  saveRDS(module_meta, file.path(time_dir, "native_module_metadata.rds"), compress = FALSE)
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
reconstruction_summary <- do.call(rbind, reconstruction_rows)
if (!identical(summary$timepoint, timepoints) || sum(summary$control_cells) != 11447L) {
  stop("Pando timepoint suite failed control-cell accounting.")
}
if (!identical(reconstruction_summary$timepoint, timepoints) ||
  any(reconstruction_summary$status != "PASS") ||
  any(reconstruction_summary$reconstructed_tf_target_edges !=
    reconstruction_summary$observed_tf_target_edges)) {
  stop("Pando timepoint suite failed module reconstruction accounting.")
}
write_tsv(summary, file.path(stage, "timepoint_run_summary.tsv"), na = "")
write_tsv(
  reconstruction_summary,
  file.path(stage, "native_module_reconstruction_suite.tsv"), na = ""
)
write_tsv(
  data.frame(
    parameter = c(
      "input_path", "input_sha256", "input_validation_sha256", "Pando_version",
      "Pando_source_commit", "network_unit", "network_pooling", "score_pooling",
      "training_cells", "heldout_cells", "native_region_filter", "exclude_exons",
      "motif_collection", "variable_feature_method", "variable_feature_nfeatures",
      "peak_to_gene_method", "upstream", "downstream", "extend", "only_tss",
      "tf_cor", "peak_cor", "fit_method", "interaction_term", "scale",
      "coefficient_p_adjust_method", "coefficient_p_adjust_scope",
      "module_p_thresh", "module_rsq_thresh", "module_nvar_thresh",
      "module_min_genes_per_module", "module_metadata_semantics", "cores",
      "cell_cap", "input_feature_cap"
    ),
    value = c(
      input_path, input_sha256, unname(sha256_files(input_validation)), as.character(packageVersion("Pando")),
      trimws(readLines(source_commit_path, warn = FALSE)[[1]]),
      "one independently fitted network per control timepoint", "none", "none",
      "all 11447 released NT controls", "all 126157 non-NT cells",
      "Pando phastConsElements20Mammals.UCSC.hg38", "true", "Pando motifs and motif2tf",
      "Seurat vst executed within each full timepoint input", "2000 (Pando variable-feature workflow)",
      "GREAT", "100000", "0", "1000000", "false", "0.1", "0", "glm", ":", "false",
      "fdr", "all fitted coefficient p-values within one independently fitted timepoint network",
      as.character(pando_p_thresh), as.character(pando_rsq_thresh),
      as.character(pando_nvar_thresh), as.character(pando_min_genes_per_module),
      "find_modules metadata after model R-squared, variable-count, and coefficient FDR gates; min_genes_per_module governs module feature lists, not metadata rows",
      as.character(cores), "none", "none"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(stage, "run_parameters.tsv"), na = ""
)
write_tsv(
  data.frame(
    field = c(
      "dataset", "author_annotated_release_cells", "training_control_cells", "heldout_perturbation_cells",
      "timepoint_networks", "complete_timepoints", "complete_replicates", "cell_cap",
      "input_feature_cap", "method_native_filtering", "native_output", "score_scope",
      "cross_timepoint_score_comparison", "causal_evaluation_scope"
    ),
    value = c(
      "GSE274113", "137604", "11447", "126157", nrow(summary),
      paste(timepoints, collapse = ";"), paste0("rep", c(1:10, 12:14, 16), collapse = ";"),
      "none", "none", "Pando phastCons plus per-timepoint variable features",
      "Pando coefficient, goodness-of-fit, independently reconstructed find_modules metadata, module features, and Network objects",
      "within_timepoint_native_endpoint_only", "prohibited",
      "replicate-paired pseudobulk effects from all heldout non-NT cells within TF-timepoint strata; author S2/S3 tables provide TF-level aggregate concordance only"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(stage, "run_validation.tsv"), na = ""
)
if (!file.rename(stage, destination)) stop("Unable to promote formal Pando result.")
cat("destination\t", destination, "\n", sep = "")
cat("completed_utc\t", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n", sep = "")
print(summary)
}



if (sys.nframe() == 0L) run_pando()
