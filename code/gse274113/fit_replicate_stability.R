#!/usr/bin/env Rscript

source("code/functions/runtime.R")
configure_task_runtime()
source("code/functions/network_io.R")
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(multiCSN)
})
data.table::setDTthreads(1L)

# Stability weights use exact evidence ties within each fit. Keep this rule
# separate from the tolerance grouping used by the layered endpoint workflow.
replicate_edge_weights <- function(beta, evidence) {
  if (length(beta) != length(evidence) || any(!is.finite(beta)) ||
      any(!is.finite(evidence)) || any(beta == 0) || any(evidence < 0)) {
    stop("Invalid replicate edge coefficient or deletion evidence.")
  }
  if (!length(beta)) return(numeric())
  ranks <- rank(-evidence, ties.method = "min")
  groups <- match(ranks, sort(unique(ranks)))
  sign(beta) * (1 - (groups - 0.5) / length(unique(groups)))
}

replicate_fraction <- function(numerator, denominator) {
  if (denominator > 0) numerator / denominator else NA_real_
}

replicate_correlation <- function(left, right) {
  if (length(left) < 3L || length(unique(left)) < 2L ||
      length(unique(right)) < 2L) return(NA_real_)
  suppressWarnings(cor(left, right, method = "spearman"))
}

replicate_node_stability <- function(full, fold, node_column) {
  full_degree <- full[, .(full_degree = .N), by = node_column]
  fold_degree <- fold[, .(fold_degree = .N), by = node_column]
  shared <- merge(
    full[, .(regulator, target)], fold[, .(regulator, target)],
    by = c("regulator", "target")
  )[, .(shared_degree = .N), by = node_column]
  result <- Reduce(function(left, right) merge(
    left, right, by = node_column, all = TRUE
  ), list(full_degree, fold_degree, shared))
  for (column in c("full_degree", "fold_degree", "shared_degree")) {
    set(result, which(is.na(result[[column]])), column, 0L)
  }
  result[, union_degree := full_degree + fold_degree - shared_degree]
  result[, jaccard := fifelse(
    union_degree > 0L, shared_degree / union_degree, NA_real_
  )]
  result
}

fit_replicate_fold <- function(
  gene_layer, gene_counts, peak_counts, regulators_declared, features,
  peaks2gene, peak_tf_gate_declared, fold_cells, full_network,
  timepoint, heldout, heldout_cells, training_replicates
) {
  fold_started <- Sys.time()
  if (length(fold_cells) < 3L || anyDuplicated(fold_cells) ||
      any(!fold_cells %in% colnames(gene_layer)) ||
      any(!features %in% rownames(gene_layer)) ||
      any(!regulators_declared %in% rownames(gene_layer))) {
    stop("Invalid replicate fit input.")
  }
  full_network <- copy(full_network)
  full_network[, full_rank := frank(-abs(weight), ties.method = "first")]
  full_indices_by_target <- split(
    seq_len(nrow(full_network)), factor(full_network$target, levels = features)
  )
  x_raw <- t(as.matrix(gene_layer[regulators_declared, fold_cells, drop = FALSE]))
  x_means <- colMeans(x_raw)
  x_centered <- sweep(x_raw, 2L, x_means, FUN = "-")
  x_scales <- sqrt(colSums(x_centered ^ 2) / (nrow(x_centered) - 1L))
  x_variable <- is.finite(x_scales) & x_scales > 0
  x_standardized <- sweep(
    x_centered[, x_variable, drop = FALSE],
    2L, x_scales[x_variable], FUN = "/"
  )
  regulators <- colnames(x_standardized)
  regulator_detected <- as.integer(Matrix::rowSums(
    gene_counts[regulators, fold_cells, drop = FALSE] != 0
  ))
  names(regulator_detected) <- regulators
  regulator_eligible <- is.finite(regulator_detected) & regulator_detected >= 20L
  gram <- crossprod(x_standardized)

  gene_data <- Matrix::t(gene_layer[features, fold_cells, drop = FALSE])
  target_detected <- as.integer(Matrix::colSums(gene_data != 0))
  names(target_detected) <- features
  target_means <- as.numeric(Matrix::colSums(gene_data)) / length(fold_cells)
  target_sumsq <- as.numeric(Matrix::colSums(gene_data ^ 2))
  target_centered_sumsq <- pmax(
    0, target_sumsq - length(fold_cells) * target_means ^ 2
  )
  target_scales <- sqrt(target_centered_sumsq / (length(fold_cells) - 1L))
  target_variable <- is.finite(target_scales) & target_scales > 0
  target_eligible <- (is.finite(target_detected) & target_detected >= 20L) & target_variable
  xty <- as.matrix(Matrix::crossprod(x_standardized, gene_data))
  xty[, target_variable] <- sweep(
    xty[, target_variable, drop = FALSE], 2L,
    target_scales[target_variable], FUN = "/"
  )
  xty[, !target_variable] <- 0
  y_ss <- ifelse(target_variable, length(fold_cells) - 1, 0)

  accessible_peak <- Matrix::rowSums(
    peak_counts[rownames(peak_tf_gate_declared), fold_cells, drop = FALSE] != 0
  ) >= 3L
  fold_peaks <- rownames(peak_tf_gate_declared)[accessible_peak]
  fold_peaks2gene <- peaks2gene[, fold_peaks, drop = FALSE]
  fold_peak_tf_gate <- peak_tf_gate_declared[
    fold_peaks, regulators, drop = FALSE
  ]

  candidate_predictors <- vector("list", length(features))
  accounting <- vector("list", length(features))
  full_candidate_present <- logical(nrow(full_network))
  full_fit_eligible <- logical(nrow(full_network))
  for (target_index in seq_along(features)) {
    target <- features[[target_index]]
    target_regions <- if (target %in% rownames(fold_peaks2gene)) {
      colnames(fold_peaks2gene)[as.logical(
        fold_peaks2gene[target, , drop = TRUE]
      )]
    } else character()
    supplied <- if (length(target_regions)) {
      gate <- fold_peak_tf_gate[target_regions, , drop = FALSE]
      which(Matrix::colSums(gate != 0) > 0)
    } else integer()
    supplied_names <- regulators[supplied]
    full_rows <- full_indices_by_target[[target_index]]
    if (length(full_rows)) {
      full_candidate_present[full_rows] <-
        full_network$regulator[full_rows] %chin% supplied_names
    }

    if (!target_eligible[[target_index]]) {
      filtered <- list(
        indices = integer(), supplied = length(supplied),
        regulator_ineligible = 0L, self = 0L, affine_alias = 0L
      )
      status <- "target_detection_or_variance_ineligible"
    } else {
      detected <- regulator_eligible[supplied]
      self <- regulators[supplied] == target
      correlation <- abs(xty[supplied, target_index]) / (length(fold_cells) - 1L)
      alias <- is.finite(correlation) & correlation >= 1 - 1e-10
      filtered <- list(
        indices = supplied[detected & !self & !alias],
        regulator_ineligible = sum(!detected), self = sum(detected & self),
        affine_alias = sum(detected & !self & alias)
      )
      status <- if (length(filtered$indices)) {
        "eligible_for_fit"
      } else "no_eligible_candidate"
    }
    candidate_predictors[[target_index]] <- filtered$indices
    filtered_names <- regulators[filtered$indices]
    if (length(full_rows)) {
      full_fit_eligible[full_rows] <-
        full_network$regulator[full_rows] %chin% filtered_names
    }
    accounting[[target_index]] <- data.frame(
      target = target,
      detected_remaining_control_cells = target_detected[[target_index]],
      target_variable = target_variable[[target_index]],
      target_eligible = target_eligible[[target_index]],
      fold_domain_candidate_tfs = length(supplied),
      removed_regulator_detection = filtered$regulator_ineligible,
      removed_self = filtered$self,
      removed_affine_alias = filtered$affine_alias,
      eligible_candidate_tfs = length(filtered$indices),
      prefit_status = status,
      stringsAsFactors = FALSE
    )
  }

  fit <- inferCSN::fit_greedy_l0_batch(
    gram = gram, xty = xty, response_ss = y_ss,
    candidates = candidate_predictors, n_obs = length(fold_cells),
    max_support_size = NULL, min_improvement = 1e-10
  )
  accounting <- rbindlist(accounting)
  accounting[, `:=`(
    solver_candidate_tfs = as.integer(fit$candidate_size),
    selected_tfs = as.integer(fit$support_size),
    rss = as.numeric(fit$rss), bic = as.numeric(fit$bic),
    solver_status = as.character(fit$status)
  )]
  if (!identical(accounting$eligible_candidate_tfs,
                 accounting$solver_candidate_tfs)) {
    stop("Candidate accounting differs from solver: ", heldout)
  }

  edge_target_index <- as.integer(fit$target_index)
  edge_predictor_index <- as.integer(fit$predictor_index)
  edges <- data.table(
    regulator = regulators[edge_predictor_index],
    target = features[edge_target_index],
    standardized_beta = as.numeric(fit$standardized_beta),
    deletion_delta_bic = as.numeric(fit$deletion_delta_bic),
    regulator_detected_control_cells = regulator_detected[edge_predictor_index],
    target_detected_control_cells = target_detected[edge_target_index]
  )
  if (!nrow(edges)) stop("Empty replicate network: ", heldout)
  edges[, weight := replicate_edge_weights(standardized_beta, deletion_delta_bic)]
  edges[, weight_magnitude__ := abs(weight)]
  setorder(edges, -weight_magnitude__, regulator, target)
  edges[, weight_magnitude__ := NULL]
  if (anyDuplicated(edges[, .(regulator, target)]) ||
      any(edges$regulator == edges$target) ||
      any(edges$regulator_detected_control_cells < 20L) ||
      any(edges$target_detected_control_cells < 20L)) {
    stop("Fold edge contract failed: ", heldout)
  }

  full_gene_set <- unique(c(full_network$regulator, full_network$target))
  remaining_detection <- as.integer(Matrix::rowSums(
    gene_counts[full_gene_set, fold_cells, drop = FALSE] != 0
  ))
  names(remaining_detection) <- full_gene_set
  regulator_variable_full <- full_network$regulator %chin% regulators
  target_variable_full <- target_variable[match(full_network$target, features)]
  rna_detection_eligible <-
    remaining_detection[full_network$regulator] >= 20L &
    remaining_detection[full_network$target] >= 20L
  rna_model_eligible <- rna_detection_eligible & regulator_variable_full &
    target_variable_full
  fold_keys <- paste(edges$regulator, edges$target, sep = "\r")
  full_keys <- paste(full_network$regulator, full_network$target, sep = "\r")
  selected_in_fold <- full_keys %chin% fold_keys
  if (any(full_fit_eligible & (!rna_model_eligible | !full_candidate_present)) ||
      any(selected_in_fold & !full_fit_eligible)) {
    stop("Full-edge decomposition contract failed: ", heldout)
  }
  decomposition <- data.table(
    regulator = full_network$regulator,
    target = full_network$target,
    full_weight = full_network$weight,
    full_rank = full_network$full_rank,
    rna_detection_eligible = rna_detection_eligible,
    regulator_variable = regulator_variable_full,
    target_variable = target_variable_full,
    rna_model_eligible = rna_model_eligible,
    multiomic_candidate_present = full_candidate_present,
    fit_eligible = full_fit_eligible,
    selected_in_fold = selected_in_fold
  )

  shared <- merge(
    full_network[, .(regulator, target, weight)],
    edges[, .(regulator, target, fold_weight = weight)],
    by = c("regulator", "target")
  )
  target_detail <- replicate_node_stability(full_network, edges, "target")
  regulator_detail <- replicate_node_stability(full_network, edges, "regulator")
  summary <- data.table(
    dataset = "GSE274113", timepoint = timepoint,
    heldout_replicate = heldout,
    heldout_cells = heldout_cells,
    training_cells = length(fold_cells),
    training_replicates = training_replicates,
    declared_regulators = length(regulators_declared),
    variable_regulators = length(regulators),
    eligible_regulators = sum(regulator_eligible),
    annotated_peaks = nrow(peak_tf_gate_declared),
    accessible_annotated_peaks = length(fold_peaks),
    fold_domain_candidate_pairs = sum(accounting$fold_domain_candidate_tfs),
    eligible_candidate_pairs = sum(accounting$eligible_candidate_tfs),
    eligible_targets = sum(target_eligible),
    full_edges = nrow(full_network), fold_edges = nrow(edges),
    shared_edges = nrow(shared),
    edge_jaccard = nrow(shared) /
      (nrow(full_network) + nrow(edges) - nrow(shared)),
    full_edge_rna_detection_eligible = sum(rna_detection_eligible),
    full_edge_rna_model_eligible = sum(rna_model_eligible),
    full_edge_candidate_present = sum(full_candidate_present),
    full_edge_candidate_and_rna_model = sum(
      full_candidate_present & rna_model_eligible
    ),
    full_edge_fit_eligible = sum(full_fit_eligible),
    full_edge_selected = sum(selected_in_fold),
    full_edge_candidate_retention_given_rna_model = replicate_fraction(
      sum(full_candidate_present & rna_model_eligible), sum(rna_model_eligible)
    ),
    full_edge_fit_eligibility_given_rna_and_candidate = replicate_fraction(
      sum(full_fit_eligible), sum(rna_model_eligible & full_candidate_present)
    ),
    full_edge_selection_given_fit_eligible = replicate_fraction(
      sum(selected_in_fold), sum(full_fit_eligible)
    ),
    shared_edge_sign_agreement = mean(
      sign(shared$weight) == sign(shared$fold_weight)
    ),
    shared_edge_weight_spearman = replicate_correlation(
      abs(shared$weight), abs(shared$fold_weight)
    ),
    median_target_jaccard = median(target_detail$jaccard, na.rm = TRUE),
    mean_target_jaccard = mean(target_detail$jaccard, na.rm = TRUE),
    median_regulator_jaccard = median(regulator_detail$jaccard, na.rm = TRUE),
    mean_regulator_jaccard = mean(regulator_detail$jaccard, na.rm = TRUE),
    elapsed_seconds = as.numeric(difftime(
      Sys.time(), fold_started, units = "secs"
    ))
  )
  list(edges = edges, accounting = accounting, decomposition = decomposition,
       target_stability = target_detail, regulator_stability = regulator_detail,
       summary = summary)
}

fit_replicate_stability <- function(
  input = "data/gse274113/objects/GSE274113_control_full/candidates/GSE274113_NT_timepoint_candidates.rds",
  validation = "data/gse274113/objects/GSE274113_control_full/candidates/candidate_validation.tsv",
  output = "results/gse274113/replicate_deletion_stability/timepoints"
) {
  # All scientific settings are fixed here and in fit_replicate_fold(): RNA
  # detection >=20 cells, peak detection >=3 cells, 100 kb upstream domains,
  # self/affine-alias exclusion, uncapped Greedy-L0+BIC and one solver thread.
  if (file.exists(output)) stop("Output already exists: ", output)
  input <- normalizePath(input, mustWork = TRUE)
  validation <- normalizePath(validation, mustWork = TRUE)
  timepoints <- c(day_7 = "day 7", day_9 = "day 9", day_11 = "day 11", day_14 = "day 14")
  expected_cells <- c(day_7 = 3304L, day_9 = 4167L, day_11 = 1903L, day_14 = 2073L)
  expected_folds <- c(day_7 = 4L, day_9 = 4L, day_11 = 3L, day_14 = 3L)
  v <- fread(validation, colClasses = "character")
  if (!all(c("field", "value") %in% names(v)) ||
      !identical(v[field == "training_control_cells", value], "11447") ||
      !identical(v[field == "filter_mode", value], "unfiltered")) {
    stop("Expected the full unfiltered 11,447-cell NT-control candidate object.")
  }
  object <- readRDS(input)
  metadata <- object[[]]
  if (ncol(object) != 11447L ||
      !all(c("Timepoint", "replicate") %in% names(metadata)) ||
      anyNA(metadata$Timepoint) || anyNA(metadata$replicate) ||
      any(!grepl("^[A-Za-z0-9_.-]+$", as.character(metadata$replicate))) ||
      !setequal(unique(as.character(metadata$Timepoint)), unname(timepoints))) {
    stop("Invalid complete-control timepoint or biological-replicate metadata.")
  }
  for (slug in names(timepoints)) {
    rows <- as.character(metadata$Timepoint) == timepoints[[slug]]
    if (sum(rows) != expected_cells[[slug]] ||
        uniqueN(as.character(metadata$replicate[rows])) != expected_folds[[slug]]) {
      stop("Timepoint cell or replicate count mismatch: ", slug)
    }
  }
  params <- multiCSN::Params(object)
  regulators_requested <- multiCSN::get_attribute(object, attribute = "tfs")
  gene_annotation <- multiCSN::get_peak_annotation(object, params$peak_assay)
  features <- intersect(unique(gene_annotation$gene_name), rownames(object[[params$rna_assay]]))
  features <- sort(features[!is.na(features) & nzchar(features)])
  regions <- multiCSN::NetworkRegions(object)
  motif2tf_all <- regions@motifs2tfs
  if (!length(features) || !length(regulators_requested)) stop("Empty target or TF catalogue.")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(".replicate_fit_", tmpdir = dirname(output))
  dir.create(stage)
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)

  for (slug in names(timepoints)) {
    timepoint <- timepoints[[slug]]
    cells <- colnames(object)[as.character(metadata$Timepoint) == timepoint]
    replicate_by_cell <- setNames(as.character(metadata[cells, "replicate"]), cells)
    replicates <- sort(unique(replicate_by_cell))
    timepoint_object <- subset(object, cells = cells)
    # Log-normalize once per complete timepoint. Fold fits then subset this
    # layer and re-standardize their own design and responses.
    timepoint_object <- Seurat::NormalizeData(timepoint_object, assay = params$rna_assay, verbose = FALSE)
    gene_layer <- multiCSN::get_layer_data(timepoint_object, assay = params$rna_assay, layer = "data")
    gene_counts <- multiCSN::get_layer_data(timepoint_object, assay = params$rna_assay, layer = "counts")
    peak_counts <- multiCSN::get_layer_data(timepoint_object, assay = params$peak_assay, layer = "counts")
    regulators_declared <- intersect(regulators_requested, rownames(gene_layer))
    peaks_requested <- multiCSN::get_attribute(object, celltypes = timepoint, attribute = "peaks")
    peak_ranges <- multiCSN::parse_peak_ranges(peaks_requested)
    names(peak_ranges) <- peaks_requested
    peaks_near_gene <- multiCSN::find_peaks_near_genes(
      peaks = peak_ranges, method = "Signac",
      genes = gene_annotation[gene_annotation$gene_name %in% features, ],
      upstream = 100000, downstream = 0, only_tss = FALSE, verbose = FALSE
    )
    peaks2gene <- multiCSN::collapse_peak_gene_domains(peaks_near_gene)
    peaks2motif <- regions@motifs@data[peaks_requested, , drop = FALSE]
    if (!identical(colnames(peaks2gene), rownames(peaks2motif))) {
      stop("Genomic-domain and motif peak order differs: ", slug)
    }
    motif_use <- intersect(colnames(peaks2motif), rownames(motif2tf_all))
    annotation_peak_use <- Matrix::colSums(peaks2gene != 0) > 0 & Matrix::rowSums(peaks2motif != 0) > 0
    peaks2gene <- peaks2gene[, annotation_peak_use, drop = FALSE]
    gate <- peaks2motif[annotation_peak_use, motif_use, drop = FALSE] %*%
      motif2tf_all[motif_use, regulators_declared, drop = FALSE]
    gate[is.na(gate)] <- 0
    gate[gate < 0] <- 0
    args <- list(gene_layer = gene_layer, gene_counts = gene_counts, peak_counts = peak_counts,
                 regulators_declared = regulators_declared, features = features,
                 peaks2gene = peaks2gene, peak_tf_gate_declared = gate, timepoint = timepoint)
    timepoint_root <- file.path(stage, slug)
    dir.create(timepoint_root)
    full_root <- file.path(timepoint_root, "detected_20")
    dir.create(full_root)
    # Fit the complete-timepoint reference with the same detection and weight
    # rules as the folds; layered endpoint tables have a different contract.
    full <- do.call(fit_replicate_fold, c(args, list(
      fold_cells = cells, full_network = data.table(regulator = character(), target = character(), weight = numeric()),
      heldout = "none", heldout_cells = 0L, training_replicates = length(replicates)
    )))
    write_tsv(full$edges, file.path(full_root, "tf_gene.tsv"))
    full_network <- full$edges
    rm(full)
    summaries <- vector("list", length(replicates))
    for (i in seq_along(replicates)) {
      heldout <- replicates[[i]]
      fold_cells <- cells[replicate_by_cell[cells] != heldout]
      fit <- do.call(fit_replicate_fold, c(args, list(
        fold_cells = fold_cells, full_network = full_network, heldout = heldout,
        heldout_cells = length(cells) - length(fold_cells), training_replicates = length(replicates) - 1L
      )))
      fold_root <- file.path(timepoint_root, "folds", heldout)
      dir.create(fold_root, recursive = TRUE)
      write_tsv(fit$edges, file.path(fold_root, "tf_gene.tsv"))
      write_tsv(fit$accounting, file.path(fold_root, "target_fit_accounting.tsv.gz"))
      write_tsv(fit$decomposition, file.path(fold_root, "full_edge_decomposition.tsv.gz"))
      write_tsv(fit$target_stability, file.path(fold_root, "target_stability.tsv.gz"))
      write_tsv(fit$regulator_stability, file.path(fold_root, "regulator_stability.tsv.gz"))
      summaries[[i]] <- fit$summary
      cat("fold_complete\t", timepoint, "\t", heldout, "\t", nrow(fit$edges), " edges\n", sep = "")
      rm(fit)
      gc(verbose = FALSE)
    }
    write_tsv(rbindlist(summaries), file.path(timepoint_root, "fold_summary.tsv"))
    rm(timepoint_object, gene_layer, gene_counts, peak_counts, gate, peaks2gene, peaks2motif, args, full_network)
    gc(verbose = FALSE)
  }
  sources <- c(input, validation, "code/gse274113/fit_replicate_stability.R",
               "code/functions/runtime.R", "code/functions/network_io.R", "code/functions/table_io.R")
  for (pkg in c("multiCSN", "inferCSN")) {
    root <- system.file(package = pkg)
    sources <- c(sources, file.path(root, c("DESCRIPTION", "NAMESPACE")),
                 list.files(file.path(root, "R"), "[.]rdb$", full.names = TRUE),
                 list.files(file.path(root, "libs"), full.names = TRUE))
  }
  sources <- normalizePath(sources, mustWork = TRUE)
  write_tsv(data.table(file = sources, bytes = file.info(sources)$size,
                      sha256 = unname(sha256_files(sources))), file.path(stage, "input_manifest.tsv"))
  write_tsv(data.table(
    field = c("dataset", "folds", "parallel_workers", "solver_threads", "rna_detection",
              "peak_detection", "upstream", "downstream", "candidate_domain", "RNA_normalization",
              "cell_sampling", "feature_sampling", "max_support_size", "min_improvement", "weight"),
    value = c("GSE274113", "14", "1", "1", "20 remaining control cells", "3 remaining control cells",
              "100000", "0", "recomputed within each replicate deletion", "complete timepoint LogNormalize; fold standardization",
              "none", "none", "uncapped", "1e-10", "signed ordinal deletion evidence; exact ties within each fit")
  ), file.path(stage, "run_contract.tsv"))
  files <- sort(list.files(stage, recursive = TRUE, full.names = TRUE))
  files <- files[!file.info(files)$isdir]
  write_tsv(data.table(file = substring(files, nchar(stage) + 2L), bytes = file.info(files)$size,
                      sha256 = unname(sha256_files(files))), file.path(stage, "output_manifest.tsv"))
  if (file.exists(output) || !file.rename(stage, output)) stop("Unable to publish complete replicate stability results.")
  cat("replicate_stability\tPASS\t", output, "\n", sep = "")
  invisible(output)
}

if (sys.nframe() == 0L) {
  options <- list(
    input = "data/gse274113/objects/GSE274113_control_full/candidates/GSE274113_NT_timepoint_candidates.rds",
    validation = "data/gse274113/objects/GSE274113_control_full/candidates/candidate_validation.tsv",
    output = "results/gse274113/replicate_deletion_stability/timepoints"
  )
  do.call(fit_replicate_stability, options)
}
