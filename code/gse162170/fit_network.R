#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(Signac)
  library(multiCSN)
})
data.table::setDTthreads(1L)
source("code/functions/network_io.R")

write_network_figure_tables <- function(fit, domain_targets, output) {
  source("code/functions/table_io.R")
  accounting <- as.data.table(fit$accounting$region_gene)
  if (anyDuplicated(accounting$target) || anyNA(accounting$target)) {
    stop("Invalid region-gene target accounting.")
  }
  annotated <- nrow(accounting)
  domain_targets <- intersect(unique(domain_targets), accounting$target)
  counts <- c(
    annotated, length(domain_targets), sum(accounting$status != "no_domain_peak"),
    sum(accounting$n_candidates > 0L),
    sum(accounting$status %chin% c("fit", "empty_support")),
    sum(accounting$n_selected > 0L)
  )
  if (annotated == 0L || any(diff(counts) > 0)) stop("Invalid target funnel.")
  funnel <- data.table(
    stage = c("Annotated targets", "Genomic domain", "Observed in ≥3 cells",
              "Variable candidates", "Fit eligible", "Selected support"),
    targets = counts, fraction_of_annotated_targets = counts / annotated,
    retention_from_previous = counts / c(annotated, head(counts, -1L))
  )
  topology <- rbindlist(lapply(c("tf_gene", "tf_region", "region_gene"), function(endpoint) {
    data.table(relation_type = endpoint,
               relation_set = c("Candidate pairs", "Selected edges"),
               relations = c(sum(fit$accounting[[endpoint]]$n_candidates),
                             nrow(fit[[endpoint]])))
  }))
  setorder(topology, relation_set, relation_type)
  status_labels <- c(
    no_domain_peak = "No eligible domain peak",
    no_variable_domain_peak = "No variable domain peak",
    constant_response = "Constant RNA response", empty_support = "Empty sparse support",
    fit = "Selected support"
  )
  barriers <- c(no_domain_peak = "candidate_absence",
                no_variable_domain_peak = "candidate_absence",
                constant_response = "response_information",
                empty_support = "sparse_selection", fit = "selected_support")
  if (any(!accounting$status %in% names(status_labels))) stop("Unknown target fit status.")
  terminal <- accounting[, .(targets = .N, candidate_pairs = sum(n_candidates),
                             selected_edges = sum(n_selected)), by = status]
  terminal[, `:=`(barrier = unname(barriers[status]),
                  fraction_of_annotated_targets = targets / annotated,
                  barrier_label = unname(status_labels[status]), fraction = targets / annotated)]
  setorder(terminal, barrier, status)
  setcolorder(terminal, c("barrier", "status", "targets", "fraction_of_annotated_targets",
                         "candidate_pairs", "selected_edges", "barrier_label", "fraction"))
  chains <- fit$triplets
  mediated <- fit$mediated_tf_gene
  values <- c(
    triplets = nrow(chains), unique_regulators = uniqueN(chains$regulator),
    unique_regions = uniqueN(chains$region), unique_targets = uniqueN(chains$target),
    positive_direction = sum(chains$chain_direction > 0),
    negative_direction = sum(chains$chain_direction < 0),
    backbone_selected = sum(chains$backbone_selected),
    backbone_selected_fraction = mean(chains$backbone_selected),
    backbone_sign_concordant = sum(chains$backbone_sign_concordant),
    backbone_sign_concordant_fraction_of_selected =
      sum(chains$backbone_sign_concordant) / sum(chains$backbone_selected),
    mediated_tf_gene = nrow(mediated),
    mediated_tied_paths = sum(mediated$strongest_path_ties > 1L),
    mediated_sign_consistent = sum(mediated$strongest_path_sign_consistent)
  )
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  write_tsv(topology, file.path(output, "network_topology.tsv"))
  write_tsv(funnel, file.path(output, "Fig4b_target_funnel.tsv"))
  write_tsv(terminal, file.path(output, "terminal_states.tsv"))
  write_tsv(data.table(measure = names(values), value = unname(values)),
            file.path(output, "chain_summary.tsv"))
}

fit_network <- function(arguments = list(
  dataset = "GSE162170",
  input = "data/gse162170/gse162170_candidates.rds",
  output = "results/gse162170/network",
  expected_cells = 8981L,
  response_chunk = 1024L
)) {
if (!"fit_layered_network" %in% getNamespaceExports("multiCSN")) {
  stop(
    "The installed multiCSN does not export fit_layered_network(). ",
    "Install the required packages using code/requirements/install_environments.sh."
  )
}

input_path <- normalizePath(arguments$input, mustWork = TRUE)
output_root <- normalizePath(arguments$output, mustWork = FALSE)
if (file.exists(output_root)) stop("Output exists and will not be overwritten: ", output_root)
dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
stage_root <- paste0(output_root, ".staging_", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"))
dir.create(stage_root, recursive = FALSE)

object <- readRDS(input_path)
cells <- multiCSN::get_attribute(object, celltypes = "aggregate", attribute = "cells")
targets_requested <- multiCSN::get_attribute(object, celltypes = "aggregate", attribute = "genes")
peaks_requested <- multiCSN::get_attribute(object, celltypes = "aggregate", attribute = "peaks")
regulators_requested <- multiCSN::get_attribute(object, attribute = "tfs")
if (length(cells) != arguments$expected_cells || anyDuplicated(cells) ||
  !length(targets_requested) || !length(peaks_requested) || !length(regulators_requested)) {
  stop("Complete candidate catalogue validation failed.")
}
# The layered multi-omic fit itself lives in the multiCSN package; this script
# only validates the catalogue, materialises the publication contract and
# promotes the staged result.
fit <- multiCSN::fit_layered_network(
  object,
  cell_group = "aggregate",
  min_detected = 0L,
  min_peak_cells = 3L,
  upstream = 100000,
  downstream = 0,
  tf_region_scope = "domain_peaks",
  exclude_response_alias = FALSE,
  renormalize = FALSE,
  response_chunk = arguments$response_chunk,
  sort_regulators = TRUE,
  checkpoint_dir = stage_root,
  verbose = TRUE
)

write_layered_network(fit, stage_root)
# Count domains before the accessibility gate using the same annotation and
# gene-body extension as the fit. Accounting supplies all subsequent gates.
params <- multiCSN::Params(object)
annotation <- Signac::Annotation(object[[params$peak_assay]])
annotation <- annotation[annotation$gene_name %in% fit$targets]
domains <- multiCSN::find_peaks_near_genes(
  peaks = Signac::StringToGRanges(peaks_requested), genes = annotation,
  method = "Signac", upstream = fit$settings$upstream,
  downstream = fit$settings$downstream, only_tss = FALSE, verbose = FALSE
)
domain_targets <- colnames(domains)[Matrix::colSums(domains != 0) > 0]
write_network_figure_tables(fit, domain_targets, stage_root)

# The component tables above are the complete result. Per-target chunk checkpoints
# only serve resume during fitting and are not part of the promoted output.
part_dirs <- list.dirs(stage_root, recursive = FALSE, full.names = TRUE)
part_dirs <- part_dirs[grepl("_parts$", part_dirs)]
if (length(part_dirs)) unlink(part_dirs, recursive = TRUE, force = TRUE)
if (!file.rename(stage_root, output_root)) stop("Unable to promote complete network result.")
cat("network_fit_complete\tPASS\t", output_root, "\n", sep = "")

}

if (sys.nframe() == 0L) fit_network()
