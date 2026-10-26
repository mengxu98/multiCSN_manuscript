#!/usr/bin/env Rscript

source("code/functions/table_io.R")

source("code/functions/runtime.R")
configure_task_runtime()

# ---- Parameters ----
# seeds and I/O paths: ablation reads multiseed outputs and writes to component_ablation.
arguments <- list(
  seeds = 1001:1010,
  input_root = "results/scmultisim/multiseed",
  output_root = "results/scmultisim/component_ablation"
)
# ---------------------------------------------------------------

repo_root <- normalizePath(getwd(), mustWork = TRUE)
component_helper <- file.path(repo_root, "code", "scmultisim", "component_helpers.R")
ranking_helper <- file.path(repo_root, "code", "functions", "ranking_metrics.R")
source(component_helper)
source(ranking_helper)

input_root <- normalizePath(arguments$input_root, mustWork = TRUE)
output_root <- normalizePath(arguments$output_root, mustWork = FALSE)
if (file.exists(output_root)) stop("Output already exists: ", output_root, call. = FALSE)
dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
stage <- paste0(output_root, ".staging_", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"))
if (file.exists(stage)) stop("Staging output already exists: ", stage, call. = FALSE)
dir.create(stage)


bic_from_rss <- function(rss, n_obs, n_terms) {
  n_obs * log(rss / n_obs) + n_terms * log(n_obs)
}

fit_full <- function(design, response, target) {
  fit <- inferCSN::fit_greedy_l0(
    x = design, y = response, max_support_size = NULL,
    min_improvement = 1e-10, verbose = FALSE
  )
  variable <- as.character(fit$coefficients$variable)
  coefficient <- as.numeric(fit$coefficients$coefficient)
  deletion_delta_bic <- as.numeric(fit$coefficients$deletion_delta_bic)
  selected <- which(is.finite(coefficient) & coefficient != 0)
  if (any(!is.finite(deletion_delta_bic[selected]))) {
    stop("Selected full-model edges lack deletion evidence.", call. = FALSE)
  }
  if (!length(selected)) {
    return(list(
      edges = data.frame(
        regulator = character(), target = character(), standardized_beta = numeric(),
        deletion_delta_bic = numeric(), stringsAsFactors = FALSE
      ),
      bic = fit$model$bic
    ))
  }
  list(
    edges = data.frame(
      regulator = variable[selected], target = target,
      standardized_beta = coefficient[selected],
      deletion_delta_bic = deletion_delta_bic[selected],
      stringsAsFactors = FALSE
    ),
    bic = fit$model$bic
  )
}

fit_forward_only <- function(design, response, target, min_improvement = 1e-10) {
  design <- scale(as.matrix(design))
  response <- as.numeric(scale(response))
  candidates <- colnames(design)
  n_obs <- nrow(design)
  selected <- character()
  current_bic <- bic_from_rss(sum(response^2), n_obs, 0L)

  repeat {
    remaining <- setdiff(candidates, selected)
    if (!length(remaining)) break
    candidate_bic <- vapply(remaining, function(candidate) {
      variables <- c(selected, candidate)
      fit <- stats::lm.fit(design[, variables, drop = FALSE], response)
      if (fit$rank < length(variables)) return(Inf)
      bic_from_rss(sum(fit$residuals^2), n_obs, length(variables))
    }, numeric(1))
    best <- remaining[[which.min(candidate_bic)]]
    if (!is.finite(candidate_bic[[best]]) ||
        candidate_bic[[best]] >= current_bic - min_improvement) break
    selected <- c(selected, best)
    current_bic <- candidate_bic[[best]]
  }

  if (!length(selected)) {
    return(list(
      edges = data.frame(regulator = character(), target = character(),
                         standardized_beta = numeric(), stringsAsFactors = FALSE),
      bic = current_bic
    ))
  }
  fit <- stats::lm.fit(design[, selected, drop = FALSE], response)
  list(
    edges = data.frame(
      regulator = selected, target = target,
      standardized_beta = as.numeric(fit$coefficients), stringsAsFactors = FALSE
    ),
    bic = current_bic
  )
}

score_edges <- function(scope, universe, edges, truth, score) {
  prediction <- edges[, c("regulator", "target"), drop = FALSE]
  prediction$score <- score
  score_factorized_scope(
    scope, universe, prediction, truth, c("regulator", "target")
  )$metrics
}

all_metrics <- list()
all_runtime <- list()
for (seed in arguments$seeds) {
  seed_started <- proc.time()[["elapsed"]]
  seed_root <- file.path(input_root, paste0("seed_", seed))
  input_path <- file.path(seed_root, "complete_simulation_input.rds")
  simulation <- readRDS(input_path)
  expression <- t(log1p(as.matrix(simulation$counts)))
  targets <- colnames(expression)
  raw_grn <- simulation$raw_grn
  regulators <- sort(unique(raw_grn$regulator))
  universe <- merge(data.frame(regulator = regulators), data.frame(target = targets))
  universe <- universe[universe$regulator != universe$target, , drop = FALSE]

  full_fits <- vector("list", length(targets))
  names(full_fits) <- targets
  full_started <- proc.time()[["elapsed"]]
  for (target in targets) {
    design <- expression[, regulators[regulators != target], drop = FALSE]
    full_fits[[target]] <- fit_full(design, expression[, target], target)
  }
  full_seconds <- proc.time()[["elapsed"]] - full_started
  forward_fits <- vector("list", length(targets))
  names(forward_fits) <- targets
  forward_started <- proc.time()[["elapsed"]]
  for (target in targets) {
    design <- expression[, regulators[regulators != target], drop = FALSE]
    forward_fits[[target]] <- fit_forward_only(design, expression[, target], target)
  }
  forward_seconds <- proc.time()[["elapsed"]] - forward_started
  full_edges <- do.call(rbind, lapply(full_fits, `[[`, "edges"))
  forward_edges <- do.call(rbind, lapply(forward_fits, `[[`, "edges"))
  if (is.null(full_edges)) full_edges <- data.frame(
    regulator = character(), target = character(), standardized_beta = numeric(),
    deletion_delta_bic = numeric(), stringsAsFactors = FALSE
  )
  if (is.null(forward_edges)) forward_edges <- data.frame(
    regulator = character(), target = character(), standardized_beta = numeric(),
    stringsAsFactors = FALSE
  )
  full_edges$component_weight <- signed_network_ordinal_weight(full_edges$standardized_beta)
  forward_edges$component_weight <- signed_network_ordinal_weight(forward_edges$standardized_beta)

  truth <- unique(raw_grn[, c("regulator", "target")])
  truth <- truth[truth$regulator != truth$target, , drop = FALSE]
  full_beta_metric <- score_edges(
    "TF-gene:static_input_GRN", universe, full_edges, truth,
    abs(full_edges$component_weight)
  )
  metric_rows <- list(
    full_beta_metric,
    score_edges("TF-gene:static_input_GRN", universe, full_edges, truth,
                abs(full_edges$deletion_delta_bic)),
    score_edges("TF-gene:static_input_GRN", universe, forward_edges, truth,
                abs(forward_edges$component_weight))
  )
  methods <- c(
    "full_add_delete_swap__beta_rank",
    "full_add_delete_swap__deletion_delta_bic_rank",
    "forward_only__beta_rank"
  )
  metrics <- do.call(rbind, metric_rows)
  metrics$method <- methods
  metrics$seed <- seed
  metrics <- metrics[, c("seed", "method", setdiff(colnames(metrics), c("seed", "method")))]
  support <- do.call(rbind, lapply(targets, function(target) {
    full_support <- full_edges$regulator[full_edges$target == target]
    forward_support <- forward_edges$regulator[forward_edges$target == target]
    data.frame(
      seed = seed, target = target,
      full_support = paste(full_support, collapse = ","),
      forward_support = paste(forward_support, collapse = ","),
      full_n_selected = length(full_support),
      forward_n_selected = length(forward_support),
      support_identical = setequal(full_support, forward_support),
      stringsAsFactors = FALSE
    )
  }))
  runtime <- data.frame(
    seed = seed, targets = length(targets),
    full_add_delete_swap_seconds = full_seconds,
    forward_only_seconds = forward_seconds,
    total_seconds = proc.time()[["elapsed"]] - seed_started,
    stringsAsFactors = FALSE
  )

  seed_output <- file.path(stage, paste0("seed_", seed))
  dir.create(seed_output)
  write_tsv(full_edges, file.path(seed_output, "full_add_delete_swap_tf_gene.tsv.gz"), writer = "table")
  write_tsv(forward_edges, file.path(seed_output, "forward_only_tf_gene.tsv.gz"), writer = "table")
  write_tsv(support, file.path(seed_output, "support_comparison.tsv"), writer = "table")
  write_tsv(metrics, file.path(seed_output, "metrics.tsv"), writer = "table")
  write_tsv(runtime, file.path(seed_output, "runtime.tsv"), writer = "table")
  all_metrics[[as.character(seed)]] <- metrics
  all_runtime[[as.character(seed)]] <- runtime
  cat("seed_completed\t", seed, "\n", sep = "")
}

metrics <- do.call(rbind, all_metrics)
runtime <- do.call(rbind, all_runtime)
summary <- do.call(rbind, lapply(split(metrics, metrics$method), function(x) {
  data.frame(
    method = x$method[[1]], seeds = nrow(x),
    auroc_median = stats::median(x$auroc),
    auroc_q25 = stats::quantile(x$auroc, 0.25, names = FALSE),
    auroc_q75 = stats::quantile(x$auroc, 0.75, names = FALSE),
    average_precision_median = stats::median(x$average_precision),
    average_precision_q25 = stats::quantile(x$average_precision, 0.25, names = FALSE),
    average_precision_q75 = stats::quantile(x$average_precision, 0.75, names = FALSE),
    n_predicted_median = stats::median(x$n_predicted),
    stringsAsFactors = FALSE
  )
}))
write_tsv(metrics, file.path(stage, "endpoint_metrics_all_seeds.tsv"), writer = "table")
write_tsv(summary, file.path(stage, "endpoint_summary.tsv"), writer = "table")
write_tsv(runtime, file.path(stage, "runtime_all_seeds.tsv"), writer = "table")
if (!file.rename(stage, output_root)) {
  stop("Could not atomically promote staging output.", call. = FALSE)
}
cat("completed\t", output_root, "\n", sep = "")
