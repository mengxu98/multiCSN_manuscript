#!/usr/bin/env Rscript

source("code/functions/table_io.R")

source("code/functions/runtime.R")
configure_task_runtime()

# ---- Parameters ----
# seeds and workers configuration (default 8 workers).
arguments <- list(
  seeds = 1001:1010,
  workers = 8L,
  output_root = "results/scmultisim/multiseed"
)
# ---------------------------------------------------------------

repo_root <- normalizePath(getwd(), mustWork = TRUE)
component_helper <- file.path(
  repo_root, "code", "scmultisim", "component_helpers.R"
)
ranking_helper <- file.path(repo_root, "code", "functions", "ranking_metrics.R")
source(component_helper)
source(ranking_helper)

output_root <- normalizePath(arguments$output_root, mustWork = FALSE)
if (file.exists(output_root)) stop("Output already exists: ", output_root)
dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
stage <- paste0(output_root, ".staging_", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"))
dir.create(stage)

suppressPackageStartupMessages(library(scMultiSim))
data("GRN_params_100", package = "scMultiSim", envir = environment())
raw_grn <- GRN_params_100
colnames(raw_grn) <- c("target", "regulator", "effect")
raw_grn$target <- paste0("g_", raw_grn$target)
raw_grn$regulator <- paste0("g_", raw_grn$regulator)

score_endpoint <- function(scope, universe, edges, truth, keys) {
  prediction <- edges[, keys, drop = FALSE]
  prediction$score <- abs(edges$component_weight)
  score_factorized_scope(scope, universe, prediction, truth, keys)
}

dynamic_truth <- function(networks, grn_parameters, released_grn) {
  mapping <- unique(rbind(
    data.frame(internal = grn_parameters$target, released = released_grn$target),
    data.frame(internal = grn_parameters$regulator, released = released_grn$regulator)
  ))
  gene_map <- stats::setNames(mapping$released, mapping$internal)
  tables <- lapply(networks, function(network) {
    index <- which(network != 0, arr.ind = TRUE)
    regulator <- unname(gene_map[paste0("gene", colnames(network)[index[, "col"]])])
    if (anyNA(regulator)) stop("Dynamic regulator ID could not be mapped.")
    data.frame(
      regulator = regulator,
      target = rownames(network)[index[, "row"]], stringsAsFactors = FALSE
    )
  })
  keys <- unlist(lapply(tables, function(x) make_key(x$regulator, x$target)))
  frequency <- sort(table(keys) / length(networks), decreasing = TRUE)
  pieces <- strsplit(names(frequency), "|", fixed = TRUE)
  data.frame(
    regulator = vapply(pieces, `[[`, character(1), 1),
    target = vapply(pieces, `[[`, character(1), 2),
    cell_frequency = as.numeric(frequency), stringsAsFactors = FALSE
  )
}

all_metrics <- list()
all_runtime <- list()
for (seed in arguments$seeds) {
  seed_started <- Sys.time()
  seed_root <- file.path(stage, paste0("seed_", seed))
  dir.create(seed_root)
  set.seed(seed)
  options <- list(
    rand.seed = seed, GRN = raw_grn, num.cells = 500L,
    tree = scMultiSim::Phyla1(), num.cifs = 50L, cif.sigma = 0.5,
    diff.cif.fraction = 0.8, do.velocity = TRUE,
    dynamic.GRN = list(
      cell.per.step = 3L, num.changing.edges = 2L,
      weight.mean = 0, weight.sd = 4
    ),
    threads = 1L, speed.up = FALSE
  )
  simulation_started <- Sys.time()
  simulation <- scMultiSim::sim_true_counts(options)
  simulation_seconds <- as.numeric(difftime(Sys.time(), simulation_started, units = "secs"))
  rna <- as.matrix(simulation$counts)
  atac <- as.matrix(simulation$atac_counts)
  region_to_gene <- as.matrix(simulation$region_to_gene)
  cell_specific_grn <- simulation$cell_specific_grn
  grn_parameters <- as.data.frame(simulation$grn_params)
  cells <- colnames(rna)
  targets <- rownames(rna)
  regulators <- sort(unique(raw_grn$regulator))
  peaks <- paste0("region_", seq_len(nrow(atac)))
  expression <- t(log1p(rna))
  accessibility_all <- t(log1p(atac))
  rownames(accessibility_all) <- cells
  colnames(accessibility_all) <- peaks
  variable <- apply(accessibility_all, 2, stats::var) > 0
  accessibility <- accessibility_all[, variable, drop = FALSE]

  fit_started <- Sys.time()
  tf_gene <- fit_component_network(
    expression[, regulators, drop = FALSE], targets, expression,
    component = "TF-gene", workers = arguments$workers, exclude_self = TRUE
  )
  region_gene <- fit_component_network(
    accessibility, targets, expression,
    component = "region-gene", workers = arguments$workers, exclude_self = FALSE
  )
  fit_seconds <- as.numeric(difftime(Sys.time(), fit_started, units = "secs"))
  colnames(tf_gene$edges)[colnames(tf_gene$edges) == "feature"] <- "regulator"
  colnames(region_gene$edges)[colnames(region_gene$edges) == "feature"] <- "peak"

  # Truth is first used after all endpoint fits have completed.
  tf_truth_static <- unique(raw_grn[, c("regulator", "target")])
  tf_truth_static <- tf_truth_static[tf_truth_static$regulator != tf_truth_static$target, ]
  tf_truth_dynamic <- dynamic_truth(cell_specific_grn, grn_parameters, raw_grn)
  tf_truth_dynamic <- tf_truth_dynamic[
    tf_truth_dynamic$regulator != tf_truth_dynamic$target, , drop = FALSE
  ]
  rg_index <- which(region_to_gene != 0, arr.ind = TRUE)
  rg_truth <- unique(data.frame(
    peak = peaks[rg_index[, "row"]], target = targets[rg_index[, "col"]],
    stringsAsFactors = FALSE
  ))
  rg_truth <- rg_truth[rg_truth$peak %in% colnames(accessibility), ]
  tf_universe <- merge(data.frame(regulator = regulators), data.frame(target = targets))
  tf_universe <- tf_universe[tf_universe$regulator != tf_universe$target, ]
  rg_universe <- merge(
    data.frame(peak = colnames(accessibility)), data.frame(target = targets)
  )

  scores <- list(
    score_endpoint(
      "TF-gene:static_input_GRN", tf_universe, tf_gene$edges,
      tf_truth_static, c("regulator", "target")
    ),
    score_endpoint(
      "TF-gene:dynamic_union", tf_universe, tf_gene$edges,
      tf_truth_dynamic[, c("regulator", "target")], c("regulator", "target")
    ),
    score_endpoint(
      "region-gene:all_variable_regions", rg_universe, region_gene$edges,
      rg_truth, c("peak", "target")
    )
  )
  metrics <- do.call(rbind, lapply(scores, `[[`, "metrics"))
  metrics$ranking_scheme <- "absolute_coefficient"
  metrics$seed <- seed
  metrics <- metrics[, c("seed", setdiff(colnames(metrics), "seed"))]
  scored_tf_static <- scores[[1]]$scored
  scored_tf_dynamic <- scores[[2]]$scored
  scored_rg <- scores[[3]]$scored
  write_tsv(tf_gene$edges, file.path(seed_root, "tf_gene.tsv.gz"), writer = "table")
  write_tsv(region_gene$edges, file.path(seed_root, "region_gene.tsv.gz"), writer = "table")
  write_tsv(rbind(tf_gene$metrics, region_gene$metrics), file.path(seed_root, "fit_metrics.tsv.gz"), writer = "table")
  write_tsv(metrics, file.path(seed_root, "endpoint_metrics.tsv"), writer = "table")
  write_tsv(scored_tf_static, file.path(seed_root, "scored_tf_gene_static.tsv.gz"), writer = "table")
  write_tsv(scored_tf_dynamic, file.path(seed_root, "scored_tf_gene_dynamic.tsv.gz"), writer = "table")
  write_tsv(scored_rg, file.path(seed_root, "scored_region_gene.tsv.gz"), writer = "table")
  write_tsv(tf_truth_dynamic, file.path(seed_root, "dynamic_tf_gene_truth.tsv.gz"), writer = "table")
  input_record <- list(
    seed = seed, options = options, raw_grn = raw_grn, counts = rna,
    atac_counts = atac, region_to_gene = region_to_gene,
    cell_specific_grn = cell_specific_grn, grn_parameters = grn_parameters
  )
  saveRDS(input_record, file.path(seed_root, "complete_simulation_input.rds"))
  runtime <- data.frame(
    seed = seed, cells = ncol(rna), genes = nrow(rna), regions = nrow(atac),
    variable_regions = ncol(accessibility), regulators = length(regulators),
    workers = arguments$workers, simulation_seconds = simulation_seconds,
    fit_seconds = fit_seconds,
    total_seconds = as.numeric(difftime(Sys.time(), seed_started, units = "secs")),
    truth_loaded_after_fit = TRUE, unsupported_endpoints = "TF-region;triplet",
    stringsAsFactors = FALSE
  )
  write_tsv(runtime, file.path(seed_root, "runtime.tsv"), writer = "table")
  all_metrics[[as.character(seed)]] <- metrics
  all_runtime[[as.character(seed)]] <- runtime
  cat("seed_completed\t", seed, "\tseconds=", runtime$total_seconds, "\n", sep = "")
}

metrics <- do.call(rbind, all_metrics)
runtime <- do.call(rbind, all_runtime)
write_tsv(metrics, file.path(stage, "endpoint_metrics_all_seeds.tsv"), writer = "table")
write_tsv(runtime, file.path(stage, "runtime_all_seeds.tsv"), writer = "table")
summary <- do.call(rbind, lapply(split(metrics, metrics$scope), function(x) {
  data.frame(
    scope = x$scope[1], seeds = nrow(x),
    n_universe_median = stats::median(x$n_universe),
    n_truth_median = stats::median(x$n_structural_truth),
    n_predicted_median = stats::median(x$n_predicted),
    auroc_median = stats::median(x$auroc),
    auroc_q25 = stats::quantile(x$auroc, 0.25, names = FALSE),
    auroc_q75 = stats::quantile(x$auroc, 0.75, names = FALSE),
    average_precision_median = stats::median(x$average_precision),
    average_precision_q25 = stats::quantile(x$average_precision, 0.25, names = FALSE),
    average_precision_q75 = stats::quantile(x$average_precision, 0.75, names = FALSE)
  )
}))
write_tsv(summary, file.path(stage, "endpoint_summary.tsv"), writer = "table")
if (!file.rename(stage, output_root)) stop("Could not publish output atomically.")
