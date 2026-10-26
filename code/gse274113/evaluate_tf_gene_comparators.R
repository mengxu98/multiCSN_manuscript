#!/usr/bin/env Rscript

source("code/functions/ranking_metrics.R")

source("code/functions/table_io.R")

source("code/functions/runtime.R")
configure_task_runtime()

suppressPackageStartupMessages({
  library(data.table)
})
data.table::setDTthreads(1L)

repo_root <- normalizePath(getwd(), mustWork = TRUE)
default_truth_root <- file.path(
  repo_root, "data", "gse274113", "heldout_pseudobulk_full"
)
default_multicsn_root <- file.path(
  repo_root, "results", "gse274113", "networks", "multiCSN"
)
default_comparator_root <- file.path(
  repo_root, "results", "gse274113", "networks"
)
truth_root <- default_truth_root
multicsn_root <- default_multicsn_root
comparator_root <- default_comparator_root
pando_root <- file.path(comparator_root, "Pando")
scenicplus_root <- file.path(comparator_root, "SCENICplus")
default_destination <- file.path(
  repo_root, "results", "gse274113", "comparison"
)
destination <- default_destination
output_parent <- dirname(destination)
timepoints <- c("day 7", "day 9", "day 11", "day 14")

required <- c(file.path(truth_root, "effect_model_summary.tsv"),
  unlist(lapply(c("day_7", "day_9", "day_11", "day_14"), function(slug) c(
    file.path(multicsn_root, slug, c("tf_gene.tsv.gz", "tf_region.tsv.gz")),
    file.path(pando_root, slug, "tf_gene.tsv.gz"),
    file.path(scenicplus_root, slug, "tf_gene.tsv.gz")))))

if (any(!file.exists(required))) {
  stop(
    "Complete heldout truth and three formal method suites are required: ",
    paste(required[!file.exists(required)], collapse = ", ")
  )
}

canonical_peak <- function(x) {
  x <- as.character(x)
  matched <- regexec("^(.+?)[-:]([0-9]+)-([0-9]+)$", x, perl = TRUE)
  pieces <- regmatches(x, matched)
  vapply(pieces, function(piece) {
    if (length(piece) != 4L) {
      return(NA_character_)
    }
    paste(piece[[2]], piece[[3]], piece[[4]], sep = ":")
  }, character(1))
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L ||
    length(unique(y[keep])) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[keep], y[keep], method = "spearman"))
}

degree_gini <- function(values) {
  values <- as.numeric(values)
  if (!length(values) || any(!is.finite(values)) || any(values < 0) || sum(values) == 0) {
    return(NA_real_)
  }
  values <- sort(values)
  count <- length(values)
  2 * sum(seq_len(count) * values) / (count * sum(values)) - (count + 1) / count
}

top_degree_share <- function(values, fraction = 0.1) {
  values <- as.numeric(values)
  if (!length(values) || any(!is.finite(values)) || any(values < 0) || sum(values) == 0) {
    return(NA_real_)
  }
  count <- max(1L, ceiling(length(values) * fraction))
  sum(sort(values, decreasing = TRUE)[seq_len(count)]) / sum(values)
}

evaluate_stratum <- function(
  truth, prediction, method, endpoint, regulator, timepoint,
  comparison_universe
) {
  truth <- as.data.table(truth)
  prediction <- as.data.table(prediction)
  if (!all(c("feature", "logFC", "FDR") %in% names(truth))) {
    stop("Truth table lacks feature, logFC or FDR for ", regulator, " at ", timepoint)
  }
  if (endpoint == "TF-region") {
    truth[, feature := canonical_peak(feature)]
    prediction[, feature := canonical_peak(feature)]
    truth <- truth[!is.na(feature)]
    prediction <- prediction[!is.na(feature)]
  }
  if (anyDuplicated(truth$feature)) stop("Duplicated truth feature in a causal stratum.")
  if (nrow(prediction)) {
    required_prediction_columns <- c(
      "feature", "relation_present", "rank_magnitude", "direction_sign",
      "score_available"
    )
    if (!all(required_prediction_columns %in% names(prediction))) {
      stop("Prediction adapter lacks explicit relation/score/direction fields for ", method)
    }
    if (any(!prediction$relation_present) ||
      any(prediction$score_available &
        (!is.finite(prediction$rank_magnitude) | prediction$rank_magnitude < 0)) ||
      any(!prediction$score_available & prediction$rank_magnitude != 0) ||
      any(is.finite(prediction$direction_sign) &
        !prediction$direction_sign %in% c(-1, 0, 1))) {
      stop("Prediction adapter has inconsistent relation/score/direction fields for ", method)
    }
    if (anyDuplicated(prediction$feature)) {
      stop(
        "Prediction adapter contains duplicate feature keys for ", method,
        ", ", regulator, ", ", timepoint, "."
      )
    }
  }
  joined <- merge(
    truth[, .(feature, logFC = as.numeric(logFC), FDR = as.numeric(FDR))],
    prediction[, .(
      feature, relation_present, rank_magnitude, direction_sign, score_available
    )],
    by = "feature", all.x = TRUE, sort = FALSE
  )
  joined[is.na(relation_present), `:=`(
    relation_present = FALSE,
    rank_magnitude = 0,
    direction_sign = NA_real_,
    score_available = FALSE
  )]
  if (endpoint == "TF-gene") joined <- joined[feature != regulator]
  joined <- joined[is.finite(logFC) & is.finite(FDR)]
  if (!nrow(joined)) stop("No evaluable features for ", method, ", ", regulator, ", ", timepoint)
  joined[, prediction_magnitude := rank_magnitude]
  joined[, direction_known := is.finite(direction_sign) & direction_sign != 0]
  joined[, predicted_perturbation_effect := fifelse(
    score_available & direction_known,
    -direction_sign * prediction_magnitude,
    0
  )]
  joined[, truth_response := FDR <= 0.05]
  joined[, predicted_nonzero := score_available & prediction_magnitude > 0]
  joined[, direction_correct := truth_response & predicted_nonzero &
    direction_known &
    sign(predicted_perturbation_effect) == sign(logFC)]
  truth_positives <- sum(joined$truth_response)
  truth_prevalence <- truth_positives / nrow(joined)
  selected_truth <- joined$truth_response & joined$predicted_nonzero
  selected_truth_direction_known <- selected_truth & joined$direction_known
  average_precision <- threshold_average_precision(
    joined$truth_response, joined$prediction_magnitude
  )
  data.frame(
    method = method,
    endpoint = endpoint,
    comparison_universe = comparison_universe,
    regulator = regulator,
    timepoint = timepoint,
    evaluable_features = nrow(joined),
    truth_responses_fdr_0_05 = truth_positives,
    truth_response_prevalence = truth_prevalence,
    native_relations = sum(joined$relation_present),
    score_available_relations = sum(joined$relation_present & joined$score_available),
    predicted_nonzero = sum(joined$predicted_nonzero),
    direction_known_nonzero = sum(joined$predicted_nonzero & joined$direction_known),
    prediction_density = mean(joined$predicted_nonzero),
    truth_response_coverage = if (truth_positives) sum(selected_truth) / truth_positives else NA_real_,
    response_average_precision = average_precision,
    response_ap_lift_over_prevalence = if (truth_prevalence > 0) {
      average_precision / truth_prevalence
    } else {
      NA_real_
    },
    response_auroc = tie_safe_auroc(joined$truth_response, joined$prediction_magnitude),
    response_precision_at_truth_k = expected_precision_at_k(
      joined$truth_response, joined$prediction_magnitude, truth_positives
    ),
    directional_precision_at_truth_k = expected_precision_at_k(
      joined$direction_correct, joined$prediction_magnitude, truth_positives
    ),
    direction_known_coverage_selected_truth = if (any(selected_truth)) {
      mean(joined$direction_known[selected_truth])
    } else {
      NA_real_
    },
    direction_accuracy_selected_truth = if (any(selected_truth_direction_known)) {
      mean(joined$direction_correct[selected_truth_direction_known])
    } else {
      NA_real_
    },
    signed_spearman_all = safe_spearman(joined$predicted_perturbation_effect, joined$logFC),
    signed_spearman_selected = safe_spearman(
      joined$predicted_perturbation_effect[joined$predicted_nonzero],
      joined$logFC[joined$predicted_nonzero]
    ),
    stringsAsFactors = FALSE
  )
}

source("code/functions/network_io.R")
load_canonical <- function(root, slug, endpoint = "tf_gene") {
  path <- file.path(root, slug, paste0(endpoint, ".tsv.gz"))
  table <- read_fit_checkpoint(path)
  feature <- if (endpoint == "tf_gene") "target" else "region"
  required <- c("regulator", feature, "rank_magnitude", "direction_sign", "score_available")
  if (!all(required %in% names(table)) || anyNA(table$score_available) ||
      any(!is.finite(table$rank_magnitude) | table$rank_magnitude < 0) ||
      any(!table$score_available & table$rank_magnitude != 0) ||
      any(is.finite(table$direction_sign) & !table$direction_sign %in% c(-1, 0, 1)) ||
      anyDuplicated(table[, c("regulator", feature), with = FALSE]))
    stop("Invalid canonical endpoint: ", path)
  table <- table[, .(regulator, feature = get(feature), relation_present = TRUE,
                     rank_magnitude, direction_sign, score_available)]
  list(table = table, path = path, input_rows = nrow(table))
}
load_multicsn <- function(slug, endpoint) {
  load_canonical(multicsn_root, slug, if (endpoint == "TF-gene") "tf_gene" else "tf_region")
}
load_pando <- function(slug) load_canonical(pando_root, slug)
load_scenicplus <- function(slug) load_canonical(scenicplus_root, slug)

effect_summary <- fread(file.path(truth_root, "effect_model_summary.tsv"))
truth_annotation_columns <- c(
  "perturbation_type", "expected_function_direction",
  "direction_validated", "validation_basis"
)
if (!all(c("target", "timepoint", "modality", truth_annotation_columns) %in%
  names(effect_summary))) {
  stop("Heldout truth lacks the prespecified perturbation-direction annotations.")
}
targets <- sort(unique(effect_summary$target))
method_design <- rbindlist(list(
  data.table(method = c("multiCSN", "Pando", "SCENIC+"), endpoint = "TF-gene"),
  data.table(method = "multiCSN", endpoint = "TF-region")
))
comparators <- c("Pando", "SCENIC+")
expected_effect_design <- CJ(
  target = targets, timepoint = timepoints, modality = c("ATAC", "RNA"),
  unique = TRUE
)
observed_effect_design <- unique(effect_summary[, .(target, timepoint, modality)])
direction_validated <- tolower(as.character(effect_summary$direction_validated)) %in%
  c("true", "1")
if (!length(targets) || anyNA(targets) || any(!nzchar(targets)) ||
  anyDuplicated(effect_summary[, .(target, timepoint, modality)]) ||
  !data.table::fsetequal(observed_effect_design, expected_effect_design) ||
  any(effect_summary$perturbation_type != "coding-sequence loss-of-function") ||
  any(effect_summary$expected_function_direction != "reduced TF function") ||
  any(!direction_validated) || any(!nzchar(effect_summary$validation_basis))) {
  stop("Heldout truth does not cover the input TF/timepoint/modality design.")
}
truth_direction_contract <- unique(effect_summary[, .(
  regulator = target, timepoint, perturbation_type,
  expected_function_direction,
  direction_validated = TRUE,
  validation_basis
)])
if (nrow(truth_direction_contract) != length(targets) * length(timepoints) ||
  anyDuplicated(truth_direction_contract[, .(regulator, timepoint)])) {
  stop("Perturbation-direction annotations disagree across truth modalities.")
}

dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
stage <- tempfile(pattern = ".staging_gse274113_comparison_", tmpdir = output_parent)
dir.create(stage)
metric_rows <- list()
adapter_rows <- list()
topology_rows <- list()
metric_cursor <- 1L
adapter_cursor <- 1L
topology_cursor <- 1L
for (timepoint in timepoints) {
  slug <- gsub(" ", "_", timepoint, fixed = TRUE)
  predictions <- list(
    multiCSN = load_multicsn(slug, "TF-gene"),
    Pando = load_pando(slug),
    `SCENIC+` = load_scenicplus(slug),
    multiCSN_region = load_multicsn(slug, "TF-region")
  )
  for (method_key in names(predictions)) {
    prediction <- predictions[[method_key]]
    method <- if (method_key == "multiCSN_region") "multiCSN" else method_key
    endpoint <- if (method_key == "multiCSN_region") "TF-region" else "TF-gene"
    adapter_rows[[adapter_cursor]] <- data.frame(
      method = method,
      endpoint = endpoint,
      timepoint = timepoint,
      input_rows = prediction$input_rows,
      input_representation = "canonical standardized network",
      collapsed_relations = nrow(prediction$table),
      score_available_edges = sum(prediction$table$score_available),
      score_unavailable_edges = sum(!prediction$table$score_available),
      nonzero_score_edges = sum(
        prediction$table$score_available & prediction$table$rank_magnitude > 0
      ),
      regulators = length(unique(prediction$table$regulator)),
      features = length(unique(prediction$table$feature)),
      stringsAsFactors = FALSE
    )
    adapter_cursor <- adapter_cursor + 1L
  }

  rna_truth_paths <- setNames(file.path(
    truth_root, "rna_effects", paste0(slug, "__", targets, ".tsv.gz")
  ), targets)
  rna_truth_by_tf <- lapply(rna_truth_paths, fread)
  for (tf in targets) {
    current_truth <- rna_truth_by_tf[[tf]]
    if (!all(c("feature", "logFC", "FDR") %in% names(current_truth)) ||
      anyDuplicated(current_truth$feature)) {
      stop("RNA truth features do not match: ", tf, " at ", timepoint)
    }
  }
  topology_eligibility <- list(
    full_truth_feature_universe = rbindlist(lapply(targets, function(tf) {
      data.table(
        regulator = tf,
        feature = setdiff(as.character(rna_truth_by_tf[[tf]]$feature), tf)
      )
    }))
  )
  if (any(vapply(topology_eligibility, function(table) {
    !nrow(table) || anyDuplicated(table[, .(regulator, feature)])
  }, logical(1)))) {
    stop("TF-specific topology eligibility construction failed at ", timepoint)
  }
  for (method_key in c("multiCSN", "Pando", "SCENIC+")) {
    network_table <- predictions[[method_key]]$table
    if (!all(c(
      "regulator", "feature", "relation_present", "rank_magnitude",
      "direction_sign", "score_available"
    ) %in%
      names(network_table))) {
      stop("Topology adapter schema changed for ", method_key)
    }
    for (universe_name in names(topology_eligibility)) {
      eligibility <- topology_eligibility[[universe_name]]
      universe_features <- sort(unique(eligibility$feature))
      current <- merge(
        network_table, eligibility,
        by = c("regulator", "feature"), all = FALSE, sort = FALSE
      )
      if (anyDuplicated(current[, .(regulator, feature)])) {
        stop("Duplicate TF-gene topology key for ", method_key, " at ", timepoint)
      }
      outdegree <- integer(length(targets))
      names(outdegree) <- targets
      if (nrow(current)) {
        observed_outdegree <- table(factor(current$regulator, levels = targets))
        outdegree[] <- as.integer(observed_outdegree)
      }
      indegree <- integer(length(universe_features))
      names(indegree) <- universe_features
      if (nrow(current)) {
        observed_indegree <- table(factor(current$feature, levels = universe_features))
        indegree[] <- as.integer(observed_indegree)
      }
      known_sign <- is.finite(current$direction_sign) & current$direction_sign != 0
      possible_edges <- nrow(eligibility)
      topology_rows[[topology_cursor]] <- data.frame(
        method = method_key,
        timepoint = timepoint,
        comparison_universe = universe_name,
        evaluation_regulators = length(targets),
        evaluation_targets = length(universe_features),
        possible_nonself_edges = possible_edges,
        native_tf_gene_edges = nrow(current),
        score_available_tf_gene_edges = sum(current$score_available),
        nonzero_score_tf_gene_edges = sum(
          current$score_available & current$rank_magnitude > 0
        ),
        native_edge_density = nrow(current) / possible_edges,
        regulators_with_outgoing_edges = sum(outdegree > 0),
        disconnected_regulators = sum(outdegree == 0),
        regulator_coverage_fraction = mean(outdegree > 0),
        targets_with_incoming_edges = sum(indegree > 0),
        disconnected_targets = sum(indegree == 0),
        target_coverage_fraction = mean(indegree > 0),
        edges_with_known_direction = sum(known_sign),
        positive_edges = sum(current$direction_sign[known_sign] > 0),
        negative_edges = sum(current$direction_sign[known_sign] < 0),
        positive_fraction_among_direction_known = if (any(known_sign)) {
          mean(current$direction_sign[known_sign] > 0)
        } else {
          NA_real_
        },
        outdegree_mean = mean(outdegree),
        outdegree_median = stats::median(outdegree),
        outdegree_max = max(outdegree),
        outdegree_gini = degree_gini(outdegree),
        top_decile_outdegree_share = top_degree_share(outdegree),
        indegree_mean = mean(indegree),
        indegree_median = stats::median(indegree),
        indegree_max = max(indegree),
        indegree_gini = degree_gini(indegree),
        top_decile_indegree_share = top_degree_share(indegree),
        topology_uses_native_relation_presence = TRUE,
        topology_eligibility_rule =
          "sum of TF-specific edgeR-tested feature pairs after self-edge removal and comparison-universe restriction",
        score_magnitude_pooling = "not_performed",
        stringsAsFactors = FALSE
      )
      topology_cursor <- topology_cursor + 1L
    }
  }

  for (tf in targets) {
    atac_truth_path <- file.path(
      truth_root, "atac_effects", paste0(slug, "__", tf, ".tsv.gz")
    )
    rna_truth <- rna_truth_by_tf[[tf]]
    atac_truth <- fread(atac_truth_path)
    for (method_key in c("multiCSN", "Pando", "SCENIC+")) {
      prediction <- predictions[[method_key]]$table[regulator == tf, .(
        feature, relation_present, rank_magnitude, direction_sign, score_available
      )]
      metric_rows[[metric_cursor]] <- evaluate_stratum(
        rna_truth, prediction, method_key, "TF-gene", tf, timepoint,
        "full_truth_feature_universe"
      )
      metric_cursor <- metric_cursor + 1L
    }
    region_prediction <- predictions$multiCSN_region$table[
      regulator == tf, .(
        feature, relation_present, rank_magnitude, direction_sign, score_available
      )
    ]
    metric_rows[[metric_cursor]] <- evaluate_stratum(
      atac_truth, region_prediction, "multiCSN", "TF-region", tf, timepoint,
      "full_truth_feature_universe"
    )
    metric_cursor <- metric_cursor + 1L
  }
}

metrics <- rbindlist(metric_rows, use.names = TRUE)
expected_stratum_design <- method_design[, CJ(regulator = targets, timepoint = timepoints),
  by = .(method, endpoint)]
expected_rows <- nrow(expected_stratum_design)
if (nrow(metrics) != expected_rows ||
  !fsetequal(metrics[, .(method, endpoint, regulator, timepoint)], expected_stratum_design) ||
  nrow(metrics[endpoint == "TF-gene"]) != length(targets) * length(timepoints) * sum(method_design$endpoint == "TF-gene") ||
  nrow(metrics[endpoint == "TF-region"]) != length(targets) * length(timepoints) * sum(method_design$endpoint == "TF-region") ||
  anyDuplicated(metrics[, .(
    method, endpoint, comparison_universe, regulator, timepoint
  )]) ||
  any(metrics$score_available_relations > metrics$native_relations) ||
  any(metrics$predicted_nonzero > metrics$score_available_relations) ||
  any(metrics$direction_known_nonzero > metrics$predicted_nonzero)) {
  stop("Method comparison did not produce all declared causal strata.")
}
unit_interval_columns <- c(
  "truth_response_prevalence", "prediction_density", "truth_response_coverage",
  "response_average_precision", "response_auroc", "response_precision_at_truth_k",
  "directional_precision_at_truth_k", "direction_known_coverage_selected_truth",
  "direction_accuracy_selected_truth"
)
for (column in unit_interval_columns) {
  values <- metrics[[column]]
  if (any(is.finite(values) & (values < 0 | values > 1))) {
    stop("Metric outside [0,1]: ", column)
  }
}
for (column in c("signed_spearman_all", "signed_spearman_selected")) {
  values <- metrics[[column]]
  if (any(is.finite(values) & (values < -1 | values > 1))) {
    stop("Correlation outside [-1,1]: ", column)
  }
}
write_tsv(metrics, file.path(stage, "stratum_metrics.tsv"))
adapter_summary <- rbindlist(adapter_rows)
if (nrow(adapter_summary) != nrow(method_design) * length(timepoints) ||
  !fsetequal(adapter_summary[, .(method, endpoint, timepoint)],
    method_design[, .(timepoint = timepoints), by = .(method, endpoint)]) ||
  anyDuplicated(adapter_summary[, .(method, endpoint, timepoint)]) ||
  any(adapter_summary$collapsed_relations > adapter_summary$input_rows) ||
  any(adapter_summary$score_available_edges +
    adapter_summary$score_unavailable_edges !=
    adapter_summary$collapsed_relations) ||
  any(adapter_summary$nonzero_score_edges > adapter_summary$score_available_edges)) {
  stop("Prediction adapter accounting failed.")
}
write_tsv(adapter_summary, file.path(stage, "prediction_adapter_summary.tsv"))
network_topology <- rbindlist(topology_rows)
topology_denominators <- network_topology[, .(
  possible_values = uniqueN(possible_nonself_edges),
  regulator_values = uniqueN(evaluation_regulators),
  target_values = uniqueN(evaluation_targets)
), by = .(timepoint, comparison_universe)]
if (nrow(network_topology) != sum(method_design$endpoint == "TF-gene") * length(timepoints) ||
  !fsetequal(network_topology[, .(method, timepoint)],
    method_design[endpoint == "TF-gene", .(timepoint = timepoints), by = method]) ||
  anyDuplicated(network_topology[, .(method, timepoint, comparison_universe)]) ||
  any(topology_denominators$possible_values != 1L) ||
  any(topology_denominators$regulator_values != 1L) ||
  any(topology_denominators$target_values != 1L) ||
  any(network_topology$possible_nonself_edges <= 0L) ||
  any(network_topology$native_tf_gene_edges > network_topology$possible_nonself_edges) ||
  any(network_topology$score_available_tf_gene_edges >
    network_topology$native_tf_gene_edges) ||
  any(network_topology$nonzero_score_tf_gene_edges >
    network_topology$score_available_tf_gene_edges) ||
  any(network_topology$edges_with_known_direction !=
    network_topology$positive_edges + network_topology$negative_edges) ||
  any(abs(network_topology$native_edge_density -
    network_topology$native_tf_gene_edges /
      network_topology$possible_nonself_edges) > 1e-12) ||
  any(abs(network_topology$outdegree_mean -
    network_topology$native_tf_gene_edges /
      network_topology$evaluation_regulators) > 1e-12) ||
  any(abs(network_topology$indegree_mean -
    network_topology$native_tf_gene_edges /
      network_topology$evaluation_targets) > 1e-12) ||
  any(network_topology$disconnected_regulators +
    network_topology$regulators_with_outgoing_edges !=
    network_topology$evaluation_regulators) ||
  any(network_topology$disconnected_targets +
    network_topology$targets_with_incoming_edges !=
    network_topology$evaluation_targets)) {
  stop("TF-gene topology accounting failed.")
}
write_tsv(network_topology, file.path(stage, "native_tf_gene_topology.tsv"))
truth_strata <- unique(metrics[, .(
  endpoint, comparison_universe, regulator, timepoint, evaluable_features,
  truth_responses_fdr_0_05, truth_response_prevalence
)])
truth_strata <- merge(
  truth_strata, truth_direction_contract,
  by = c("regulator", "timepoint"), all.x = TRUE, sort = FALSE
)
setcolorder(truth_strata, c(
  "endpoint", "comparison_universe", "regulator", "timepoint",
  "perturbation_type", "expected_function_direction", "direction_validated",
  "validation_basis", "evaluable_features", "truth_responses_fdr_0_05",
  "truth_response_prevalence"
))
if (nrow(truth_strata) != length(targets) * length(timepoints) * length(unique(expected_effect_design$modality)) ||
  any(!truth_strata$direction_validated) ||
  any(truth_strata$perturbation_type != "coding-sequence loss-of-function") ||
  any(truth_strata$expected_function_direction != "reduced TF function") ||
  any(!nzchar(truth_strata$validation_basis))) {
  stop("Causal truth universe differs between method adapters.")
}
write_tsv(truth_strata, file.path(stage, "truth_stratum_summary.tsv"))

metric_columns <- c(
  "prediction_density", "truth_response_coverage", "response_average_precision",
  "response_ap_lift_over_prevalence", "response_auroc",
  "response_precision_at_truth_k", "directional_precision_at_truth_k",
  "direction_known_coverage_selected_truth", "direction_accuracy_selected_truth",
  "signed_spearman_all", "signed_spearman_selected"
)
set.seed(20260803)
bootstrap_replicates <- 5000L
macro_rows <- list()
cursor <- 1L
for (current_method in unique(metrics$method)) {
  for (current_endpoint in unique(metrics[method == current_method]$endpoint)) {
    current_universes <- unique(metrics[
      method == current_method & endpoint == current_endpoint,
      comparison_universe
    ])
    for (current_universe in current_universes) {
      current <- metrics[
        method == current_method & endpoint == current_endpoint &
          comparison_universe == current_universe
      ]
      tf_clusters <- sort(unique(current$regulator))
      for (metric in metric_columns) {
        values <- current[[metric]]
        finite <- is.finite(values)
        cluster_stats <- current[, .(
          finite_sum = sum(get(metric)[is.finite(get(metric))]),
          finite_n = sum(is.finite(get(metric)))
        ), by = regulator][match(tf_clusters, regulator)]
        finite_clusters <- which(cluster_stats$finite_n > 0L)
        if (length(finite_clusters)) {
          sampled <- matrix(
            sample(
              finite_clusters, length(finite_clusters) * bootstrap_replicates,
              replace = TRUE
            ),
            nrow = length(finite_clusters)
          )
          boot_sum <- colSums(matrix(
            cluster_stats$finite_sum[sampled],
            nrow = length(finite_clusters)
          ))
          boot_n <- colSums(matrix(
            cluster_stats$finite_n[sampled],
            nrow = length(finite_clusters)
          ))
          boot <- boot_sum / boot_n
          ci_lower <- unname(stats::quantile(boot, 0.025, na.rm = TRUE))
          ci_upper <- unname(stats::quantile(boot, 0.975, na.rm = TRUE))
        } else {
          ci_lower <- NA_real_
          ci_upper <- NA_real_
        }
        macro_rows[[cursor]] <- data.frame(
          method = current_method, endpoint = current_endpoint,
          comparison_universe = current_universe, metric = metric,
          strata_total = nrow(current), strata_finite = sum(finite),
          finite_strata_fraction = sum(finite) / nrow(current),
          regulator_clusters_total = length(tf_clusters),
          regulator_clusters_finite = length(finite_clusters),
          macro_mean = if (any(finite)) mean(values[finite]) else NA_real_,
          macro_median = if (any(finite)) stats::median(values[finite]) else NA_real_,
          macro_sd = if (sum(finite) > 1L) stats::sd(values[finite]) else NA_real_,
          cluster_bootstrap_ci_lower = ci_lower,
          cluster_bootstrap_ci_upper = ci_upper,
          bootstrap_replicates = if (length(finite_clusters)) bootstrap_replicates else 0L,
          bootstrap_unit = "perturbed TF with all four timepoints",
          comparison_role =
            "descriptive within-method summary; not primary between-method evidence",
          eligible_as_primary_between_method_evidence = FALSE,
          undefined_reason = if (any(finite)) {
            NA_character_
          } else {
            "metric undefined in every TF-timepoint stratum"
          },
          stringsAsFactors = FALSE
        )
        cursor <- cursor + 1L
      }
    }
  }
}
macro_metrics <- rbindlist(macro_rows)
if (nrow(macro_metrics) != nrow(method_design) * length(metric_columns) ||
  anyDuplicated(macro_metrics[, .(
    method, endpoint, comparison_universe, metric
  )]) ||
  any(macro_metrics$strata_finite > macro_metrics$strata_total) ||
  any(abs(macro_metrics$finite_strata_fraction -
    macro_metrics$strata_finite / macro_metrics$strata_total) > 1e-12) ||
  any(macro_metrics$eligible_as_primary_between_method_evidence) ||
  any(macro_metrics$regulator_clusters_finite >
    macro_metrics$regulator_clusters_total)) {
  stop("Macro metric table is incomplete or internally inconsistent.")
}
write_tsv(macro_metrics, file.path(stage, "macro_metrics.tsv"))

timepoint_metrics <- metrics[
  , lapply(.SD, finite_mean),
  by = .(method, endpoint, comparison_universe, timepoint), .SDcols = metric_columns
]
tf_metrics <- metrics[
  , lapply(.SD, finite_mean),
  by = .(method, endpoint, comparison_universe, regulator), .SDcols = metric_columns
]
if (nrow(timepoint_metrics) != nrow(method_design) * length(timepoints) ||
  nrow(tf_metrics) != nrow(method_design) * length(targets) ||
  anyDuplicated(timepoint_metrics[, .(
    method, endpoint, comparison_universe, timepoint
  )]) ||
  anyDuplicated(tf_metrics[, .(
    method, endpoint, comparison_universe, regulator
  )])) {
  stop("Timepoint or regulator macro table is incomplete.")
}
write_tsv(timepoint_metrics, file.path(stage, "timepoint_macro_metrics.tsv"))
write_tsv(tf_metrics, file.path(stage, "regulator_macro_metrics.tsv"))

paired_metrics <- c(
  "truth_response_coverage", "response_average_precision", "response_auroc",
  "response_ap_lift_over_prevalence", "response_precision_at_truth_k",
  "directional_precision_at_truth_k",
  "direction_accuracy_selected_truth", "signed_spearman_all"
)
comparison_rows <- list()
cursor <- 1L
for (comparator in comparators) {
  for (current_universe in "full_truth_feature_universe") {
    base <- metrics[
      method == "multiCSN" & endpoint == "TF-gene" &
        comparison_universe == current_universe
    ]
    other <- metrics[
      method == comparator & endpoint == "TF-gene" &
        comparison_universe == current_universe
    ]
    paired <- merge(
      base, other,
      by = c("regulator", "timepoint", "endpoint", "comparison_universe"),
      suffixes = c("_multicsn", "_comparator"), sort = TRUE
    )
    if (nrow(paired) != length(targets) * length(timepoints)) {
      stop("Incomplete paired comparison for ", comparator, " in ", current_universe)
    }
    tf_clusters <- sort(unique(paired$regulator))
    for (metric in paired_metrics) {
      left <- paired[[paste0(metric, "_multicsn")]]
      right <- paired[[paste0(metric, "_comparator")]]
      difference <- right - left
      finite <- is.finite(difference)
      if (!any(finite)) {
        comparison_rows[[cursor]] <- data.frame(
          comparator = comparator,
          reference = "multiCSN",
          endpoint = "TF-gene",
          comparison_universe = current_universe,
          metric = metric,
          difference_definition = "comparator minus multiCSN",
          paired_strata_total = nrow(paired),
          paired_strata = 0L,
          paired_strata_fraction = 0,
          regulator_clusters_total = length(tf_clusters),
          regulator_clusters_finite = 0L,
          mean_difference = NA_real_,
          median_difference = NA_real_,
          comparator_wins = 0L,
          exact_ties = 0L,
          comparator_losses = 0L,
          cluster_bootstrap_ci_lower = NA_real_,
          cluster_bootstrap_ci_upper = NA_real_,
          bootstrap_replicates = 0L,
          comparison_role = "primary between-method evidence on jointly finite same strata",
          same_stratum_requirement = TRUE,
          stringsAsFactors = FALSE
        )
        cursor <- cursor + 1L
        next
      }
      cluster_stats <- data.table(
        regulator = paired$regulator, difference = difference
      )[is.finite(difference), .(finite_sum = sum(difference), finite_n = .N), by = regulator]
      finite_tf_clusters <- cluster_stats$regulator
      sampled <- matrix(
        sample.int(
          length(finite_tf_clusters),
          length(finite_tf_clusters) * bootstrap_replicates,
          replace = TRUE
        ),
        nrow = length(finite_tf_clusters)
      )
      boot_sum <- colSums(matrix(
        cluster_stats$finite_sum[sampled],
        nrow = length(finite_tf_clusters)
      ))
      boot_n <- colSums(matrix(
        cluster_stats$finite_n[sampled],
        nrow = length(finite_tf_clusters)
      ))
      boot <- boot_sum / boot_n
      comparison_rows[[cursor]] <- data.frame(
        comparator = comparator,
        reference = "multiCSN",
        endpoint = "TF-gene",
        comparison_universe = current_universe,
        metric = metric,
        difference_definition = "comparator minus multiCSN",
        paired_strata_total = nrow(paired),
        paired_strata = sum(finite),
        paired_strata_fraction = sum(finite) / nrow(paired),
        regulator_clusters_total = length(tf_clusters),
        regulator_clusters_finite = length(finite_tf_clusters),
        mean_difference = mean(difference[finite]),
        median_difference = stats::median(difference[finite]),
        comparator_wins = sum(difference[finite] > 0),
        exact_ties = sum(difference[finite] == 0),
        comparator_losses = sum(difference[finite] < 0),
        cluster_bootstrap_ci_lower = unname(stats::quantile(boot, 0.025, na.rm = TRUE)),
        cluster_bootstrap_ci_upper = unname(stats::quantile(boot, 0.975, na.rm = TRUE)),
        bootstrap_replicates = bootstrap_replicates,
        comparison_role = "primary between-method evidence on jointly finite same strata",
        same_stratum_requirement = TRUE,
        stringsAsFactors = FALSE
      )
      cursor <- cursor + 1L
    }
  }
}
paired_comparisons <- rbindlist(comparison_rows)
if (nrow(paired_comparisons) != length(comparators) * length(paired_metrics) ||
  anyDuplicated(paired_comparisons[, .(
    comparator, comparison_universe, metric
  )]) ||
  any(paired_comparisons$regulator_clusters_finite >
    paired_comparisons$regulator_clusters_total) ||
  any(paired_comparisons$paired_strata > paired_comparisons$paired_strata_total) ||
  any(abs(paired_comparisons$paired_strata_fraction -
    paired_comparisons$paired_strata /
      paired_comparisons$paired_strata_total) > 1e-12) ||
  any(!paired_comparisons$same_stratum_requirement) ||
  any(paired_comparisons$comparator_wins + paired_comparisons$exact_ties +
    paired_comparisons$comparator_losses != paired_comparisons$paired_strata)) {
  stop("Paired method-comparison table is incomplete or inconsistent.")
}
write_tsv(paired_comparisons, file.path(stage, "paired_method_differences.tsv"))

publish_result_files(stage, destination)
print(macro_metrics)
