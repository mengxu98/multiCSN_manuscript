empirical_upper_tail_midrank <- function(values) {
  values <- abs(suppressWarnings(as.numeric(values)))
  if (!length(values) || any(!is.finite(values)) || any(values <= 0)) {
    stop("Tail ranks require finite, non-zero magnitudes.", call. = FALSE)
  }
  (rank(-values, ties.method = "average") - 0.5) / length(values)
}

make_key <- function(...) paste(..., sep = "|")

fit_component_network <- function(
    design,
    targets,
    response,
    component,
    workers = 1L,
    exclude_self = FALSE) {
  design <- as.matrix(design)
  response <- as.matrix(response)
  if (!identical(rownames(design), rownames(response))) {
    stop("Design and response cell orders differ.", call. = FALSE)
  }
  if (!length(colnames(design)) || !length(colnames(response))) {
    stop("Component matrices require feature names.", call. = FALSE)
  }

  fit_one_target <- function(target) {
    started <- proc.time()[["elapsed"]]
    keep <- rep(TRUE, ncol(design))
    if (exclude_self) {
      keep <- colnames(design) != target
    }
    target_design <- design[, keep, drop = FALSE]
    fit <- inferCSN::fit_greedy_l0(
      x = target_design,
      y = response[, target],
      max_support_size = NULL,
      min_improvement = 1e-10,
      verbose = FALSE
    )
    coefficient <- suppressWarnings(as.numeric(fit$coefficients$coefficient))
    selected <- which(is.finite(coefficient) & coefficient != 0)
    edge <- if (length(selected)) {
      data.frame(
        feature = colnames(target_design)[selected],
        target = target,
        standardized_beta = coefficient[selected],
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        feature = character(), target = character(),
        standardized_beta = numeric(), stringsAsFactors = FALSE
      )
    }
    metric <- data.frame(
      component = component,
      target = target,
      n_cells = nrow(target_design),
      n_candidates = ncol(target_design),
      n_selected = nrow(edge),
      r_squared = fit$metrics$r_squared,
      rss = fit$model$rss,
      bic = fit$model$bic,
      seconds = proc.time()[["elapsed"]] - started,
      stringsAsFactors = FALSE
    )
    list(edge = edge, metric = metric)
  }

  workers <- min(as.integer(workers), length(targets))
  if (.Platform$OS.type == "unix" && workers > 1L) {
    fits <- parallel::mclapply(
      targets,
      fit_one_target,
      mc.cores = workers,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    fits <- lapply(targets, fit_one_target)
  }
  edges <- do.call(rbind, lapply(fits, `[[`, "edge"))
  metrics <- do.call(rbind, lapply(fits, `[[`, "metric"))
  if (is.null(edges)) {
    edges <- data.frame(
      feature = character(), target = character(),
      standardized_beta = numeric(), stringsAsFactors = FALSE
    )
  }
  rownames(edges) <- NULL
  rownames(metrics) <- NULL
  if (nrow(edges)) {
    edges$tail_rank <- empirical_upper_tail_midrank(edges$standardized_beta)
    edges$component_weight <- signed_network_ordinal_weight(edges$standardized_beta)
  } else {
    edges$tail_rank <- numeric()
    edges$component_weight <- numeric()
  }
  list(edges = edges, metrics = metrics)
}

score_factorized_scope <- function(
    scope,
    universe,
    prediction,
    truth,
    columns,
    structural_truth_n = nrow(truth)) {
  universe$key <- do.call(make_key, universe[columns])
  prediction$key <- if (nrow(prediction)) {
    do.call(make_key, prediction[columns])
  } else {
    character()
  }
  truth$key <- if (nrow(truth)) do.call(make_key, truth[columns]) else character()
  if (anyDuplicated(universe$key) || anyDuplicated(prediction$key) ||
      !all(unique(truth$key) %in% universe$key)) {
    stop("Invalid identifiers in scoring scope: ", scope, call. = FALSE)
  }
  score_by_key <- if (nrow(prediction)) {
    stats::setNames(prediction$score, prediction$key)
  } else {
    numeric()
  }
  universe$score <- unname(score_by_key[universe$key])
  universe$score[is.na(universe$score)] <- 0
  universe$truth <- universe$key %in% unique(truth$key)
  ranked <- ranking_metrics_tie_aware(
    scores = universe$score,
    labels = universe$truth,
    top_k = c(25L, 50L, 100L, 152L, 250L, 500L, 1000L, 2000L)
  )
  metrics <- data.frame(
    scope = scope,
    n_universe = nrow(universe),
    n_structural_truth = structural_truth_n,
    n_evaluable_truth = sum(universe$truth),
    observable_truth_fraction = if (structural_truth_n > 0) {
      sum(universe$truth) / structural_truth_n
    } else {
      NA_real_
    },
    prevalence = mean(universe$truth),
    n_predicted = sum(universe$score > 0),
    auroc = ranked$auroc,
    average_precision = ranked$average_precision,
    trapezoidal_auprc = ranked$trapezoidal_auprc,
    tie_policy = "score_block_average_precision",
    stringsAsFactors = FALSE
  )
  topk <- cbind(
    data.frame(scope = scope, stringsAsFactors = FALSE),
    setNames(
      ranked$topk,
      c("top_k", "n_selected", "cutoff_score", "precision", "recall", "f1")
    )
  )
  scored <- universe[
    ranked$order,
    c(columns, "score", "truth"),
    drop = FALSE
  ]
  list(metrics = metrics, topk = topk, scored = scored)
}
