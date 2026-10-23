#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
})
summary_threads <- 4L  # Worker threads
if (!is.finite(summary_threads) || summary_threads < 1L) summary_threads <- 4L
setDTthreads(summary_threads)


# ---- Parameters ----
arguments <- list(
  "results/gse194122/networks",
  "runs/gse194122/donor_deletion_runs/results",
  "results/gse194122/sampling/donor_deletion_plan.tsv",
  "results/gse194122/trajectory",
  "results/gse194122/donor_stability"
)
# ---------------------------------------------------------------
full_root <- arguments[[1L]]
jackknife_root <- arguments[[2L]]
task_path <- arguments[[3L]]
expression_root <- arguments[[4L]]
output <- arguments[[5L]]
full_root <- normalizePath(full_root, mustWork = TRUE)
jackknife_root <- normalizePath(jackknife_root, mustWork = TRUE)
task_path <- normalizePath(task_path, mustWork = TRUE)
expression_root <- normalizePath(expression_root, mustWork = TRUE)
if (file.exists(output)) stop("Output already exists: ", output, call. = FALSE)
source("code/functions/network_contract.R")
stage <- paste0(output, ".staging")
if (file.exists(stage)) stop("Staging path already exists: ", stage, call. = FALSE)
dir.create(stage, recursive = TRUE)

states <- data.table(
  state_order = 1:5,
  state = c("HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast"),
  state_slug = c("01_HSC", "02_MK_E_prog", "03_Proerythroblast", "04_Erythroblast", "05_Normoblast")
)
tasks <- fread(task_path)
if (nrow(tasks) != 49L || anyDuplicated(tasks$task_slug)) stop("Task plan changed.")
missing_tasks <- tasks[!dir.exists(file.path(jackknife_root, task_slug)), task_slug]
if (length(missing_tasks)) stop("Missing completed folds: ", paste(missing_tasks, collapse = ", "))
invisible(lapply(file.path(jackknife_root, tasks$task_slug), validate_network_result))
invisible(lapply(file.path(full_root, states$state_slug), validate_network_result))

endpoint_spec <- list(
  `TF-gene` = list(file = "tf_gene.tsv", keys = c("regulator", "target"), direction = "standardized_beta"),
  `factorized-chain` = list(file = "triplets.tsv", keys = c("regulator", "region", "target"), direction = "chain_direction")
)
read_endpoint <- function(path, endpoint) {
  multiCSN::read_network_endpoints(
    path, endpoints = endpoint, direction = TRUE,
    weight_column = endpoint_spec[[endpoint]]$direction
  )[[endpoint]]
}

support_tables <- list()
stable_tables <- list()
fold_tables <- list()
for (endpoint in names(endpoint_spec)) {
  endpoint_folds <- list()
  endpoint_support <- list()
  endpoint_stable <- list()
  for (state_index in seq_len(nrow(states))) {
    state_row <- states[state_index]
    state_tasks <- tasks[state == state_row$state]
    fold_rows <- vector("list", nrow(state_tasks))
    for (fold_index in seq_len(nrow(state_tasks))) {
      task <- state_tasks[fold_index]
      edge_table <- read_endpoint(file.path(jackknife_root, task$task_slug), endpoint)
      edge_table[, `:=`(
        state_order = state_row$state_order,
        state = state_row$state,
        heldout_donor = task$heldout_donor,
        task_slug = task$task_slug
      )]
      fold_rows[[fold_index]] <- edge_table
    }
    folded <- rbindlist(fold_rows, use.names = TRUE, fill = TRUE)
    endpoint_folds[[state_index]] <- folded
    keys <- endpoint_spec[[endpoint]]$keys
    support <- folded[, .(signed_folds = uniqueN(task_slug)), by = c(keys, "direction")]
    support[, `:=`(
      endpoint = endpoint,
      state_order = state_row$state_order,
      state = state_row$state,
      total_folds = nrow(state_tasks),
      signed_support = signed_folds / nrow(state_tasks)
    )]
    full <- read_endpoint(file.path(full_root, state_row$state_slug), endpoint)
    full[, full_selected := TRUE]
    support <- merge(support, full, by = c(keys, "direction"), all.x = TRUE, sort = FALSE)
    support[is.na(full_selected), full_selected := FALSE]
    support[, stable_consensus := full_selected & signed_support >= 0.8]
    endpoint_support[[state_index]] <- support
    endpoint_stable[[state_index]] <- support[stable_consensus == TRUE]
  }
  fold_tables[[endpoint]] <- rbindlist(endpoint_folds, fill = TRUE)
  support_tables[[endpoint]] <- rbindlist(endpoint_support, fill = TRUE)
  stable_tables[[endpoint]] <- rbindlist(endpoint_stable, fill = TRUE)
}
support_all <- rbindlist(support_tables, fill = TRUE)
stable_all <- rbindlist(stable_tables, fill = TRUE)
stability_thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
threshold_summary <- rbindlist(lapply(stability_thresholds, function(threshold) {
  support_all[full_selected == TRUE, .(
    full_selected_edges = .N,
    consensus_edges = sum(signed_support >= threshold)
  ), by = .(endpoint, state_order, state)][, `:=`(
    signed_support_threshold = threshold,
    consensus_fraction_of_full = consensus_edges / full_selected_edges
  )]
}))
setorder(threshold_summary, endpoint, state_order, signed_support_threshold)

classify_transition <- function(previous, next_state, keys, endpoint, transition_order) {
  merged <- multiCSN::classify_edge_transitions(previous, next_state, keys)
  merged[, `:=`(
    endpoint = endpoint,
    transition_order = transition_order,
    previous_state = states$state[[transition_order]],
    next_state = states$state[[transition_order + 1L]],
    transition = paste0(states$state[[transition_order]], " -> ", states$state[[transition_order + 1L]])
  )]
  merged[]
}

transition_rows <- list()
transition_index <- 0L
for (endpoint in names(endpoint_spec)) {
  keys <- endpoint_spec[[endpoint]]$keys
  stable <- stable_tables[[endpoint]]
  for (transition_order in 1:4) {
    transition_index <- transition_index + 1L
    transition_rows[[transition_index]] <- classify_transition(
      stable[state_order == transition_order],
      stable[state_order == transition_order + 1L],
      keys, endpoint, transition_order
    )
  }
}
transition_edges <- rbindlist(transition_rows, fill = TRUE)
transition_summary <- transition_edges[, .(edges = .N), by = .(
  endpoint, transition_order, previous_state, next_state, transition, topology_class
)]
transition_summary[, union_edges := sum(edges), by = .(endpoint, transition_order)]
transition_summary[, fraction := edges / union_edges]
setorder(transition_summary, endpoint, transition_order, topology_class)

regulator_summary <- transition_edges[, .(edges = .N), by = .(
  endpoint, transition_order, previous_state, next_state, transition,
  regulator, topology_class
)]
regulator_wide <- dcast(
  regulator_summary,
  endpoint + transition_order + previous_state + next_state + transition + regulator ~ topology_class,
  value.var = "edges", fill = 0
)
for (column in c("gained", "lost", "shared_sign_stable", "shared_sign_flip")) {
  if (!column %in% names(regulator_wide)) regulator_wide[, (column) := 0L]
}
regulator_wide[, `:=`(
  turnover = gained + lost + shared_sign_flip,
  net_edge_change = gained - lost,
  union_edges = gained + lost + shared_sign_stable + shared_sign_flip
)]
regulator_wide[, turnover_fraction := fifelse(union_edges > 0, turnover / union_edges, NA_real_)]
setorder(regulator_wide, endpoint, transition_order, -turnover, regulator)

tf_concordance <- merge(
  regulator_wide[endpoint == "TF-gene", .(
    transition_order, transition, regulator,
    tf_gene_gained = gained, tf_gene_lost = lost,
    tf_gene_turnover = turnover, tf_gene_turnover_fraction = turnover_fraction
  )],
  regulator_wide[endpoint == "factorized-chain", .(
    transition_order, transition, regulator,
    chain_gained = gained, chain_lost = lost,
    chain_turnover = turnover, chain_turnover_fraction = turnover_fraction
  )],
  by = c("transition_order", "transition", "regulator"), all = TRUE, sort = FALSE
)
for (column in setdiff(names(tf_concordance), c("transition_order", "transition", "regulator"))) {
  set(tf_concordance, which(is.na(tf_concordance[[column]])), column, 0)
}
tf_concordance[, concordant_turnover_score :=
  sqrt(tf_gene_turnover_fraction * chain_turnover_fraction)]
setorder(tf_concordance, transition_order, -concordant_turnover_score, -tf_gene_turnover, -chain_turnover, regulator)

activity <- stable_all[, .(
  stable_supports = .N,
  activating_supports = sum(direction > 0),
  repressing_supports = sum(direction < 0),
  net_direction = sum(direction),
  median_signed_support = median(signed_support),
  targets = uniqueN(target)
), by = .(endpoint, state_order, state, regulator)]
activity[, normalized_outdegree := stable_supports / sum(stable_supports), by = .(endpoint, state_order)]
setorder(activity, endpoint, regulator, state_order)

expression <- fread(file.path(expression_root, "donor_state_expression.tsv.gz"))
expression_groups <- fread(file.path(expression_root, "donor_state_group_counts.tsv"))
if (nrow(expression_groups) != 49L || uniqueN(expression_groups$donor_id) != 10L ||
  !all(c(
    "state", "donor_id", "feature", "feature_role",
    "mean_normalized_expression", "detection_fraction"
  ) %in% names(expression))) {
  stop("Donor-state expression schema changed.")
}
expression[, state_order := match(state, states$state)]
tf_expression <- expression[feature_role %chin% c("TF", "TF_and_target")]
tf_state_expression <- tf_expression[, .(
  donor_groups = .N,
  median_donor_expression = median(mean_normalized_expression),
  q25_donor_expression = quantile(mean_normalized_expression, 0.25),
  q75_donor_expression = quantile(mean_normalized_expression, 0.75),
  median_detection_fraction = median(detection_fraction)
), by = .(state_order, state, regulator = feature)]
activity_expression <- merge(
  activity, tf_state_expression,
  by = c("state_order", "state", "regulator"), all.x = TRUE, sort = FALSE
)
setorder(activity_expression, endpoint, regulator, state_order)

expression_transition_rows <- list()
for (tr in 1:4) {
  previous <- expression[state_order == tr, .(
    donor_id, feature,
    previous_expression = mean_normalized_expression,
    previous_detection = detection_fraction
  )]
  next_state <- expression[state_order == tr + 1L, .(
    donor_id, feature,
    next_expression = mean_normalized_expression,
    next_detection = detection_fraction
  )]
  paired <- merge(previous, next_state, by = c("donor_id", "feature"), sort = FALSE)
  paired[, `:=`(
    expression_change = next_expression - previous_expression,
    detection_change = next_detection - previous_detection
  )]
  expression_transition_rows[[tr]] <- paired[, .(
    donor_pairs = .N,
    median_previous_expression = median(previous_expression),
    median_next_expression = median(next_expression),
    median_expression_change = median(expression_change),
    q25_expression_change = quantile(expression_change, 0.25),
    q75_expression_change = quantile(expression_change, 0.75),
    positive_change_fraction = mean(expression_change > 0),
    negative_change_fraction = mean(expression_change < 0),
    median_detection_change = median(detection_change)
  ), by = feature][, `:=`(
    transition_order = tr,
    transition = paste0(
      states$state[[tr]], " -> ", states$state[[tr + 1L]]
    )
  )]
}
expression_transitions <- rbindlist(expression_transition_rows)
tf_expression_transitions <- expression_transitions[feature %in% unique(tf_expression$feature)]
tf_candidates <- merge(
  tf_concordance,
  tf_expression_transitions,
  by.x = c("transition_order", "transition", "regulator"),
  by.y = c("transition_order", "transition", "feature"),
  all.x = TRUE, sort = FALSE
)
tf_candidates[, `:=`(
  transition_turnover_q75 = quantile(concordant_turnover_score, 0.75, na.rm = TRUE),
  transition_abs_expression_q25 = quantile(abs(median_expression_change), 0.25, na.rm = TRUE)
), by = transition_order]
tf_candidates[, expression_network_decoupled :=
  donor_pairs >= 8L &
    concordant_turnover_score >= transition_turnover_q75 &
    abs(median_expression_change) <= transition_abs_expression_q25]
setorder(
  tf_candidates, transition_order, -concordant_turnover_score,
  -expression_network_decoupled, -tf_gene_turnover, -chain_turnover, regulator
)

rewired_edges <- transition_edges[
  topology_class %chin% c("gained", "lost", "shared_sign_flip")
]
target_modules <- rbindlist(list(
  rewired_edges[endpoint == "TF-gene", .(
    supports = .N, regions = NA_integer_
  ), by = .(endpoint, transition_order, transition, regulator, target, topology_class)],
  rewired_edges[endpoint == "factorized-chain", .(
    supports = .N, regions = uniqueN(region[!is.na(region)])
  ), by = .(endpoint, transition_order, transition, regulator, target, topology_class)]
), use.names = TRUE)
setorder(target_modules, endpoint, transition_order, regulator, topology_class, -supports, target)
target_modules_expression <- merge(
  target_modules,
  expression_transitions,
  by.x = c("transition_order", "transition", "target"),
  by.y = c("transition_order", "transition", "feature"),
  all.x = TRUE, sort = FALSE
)
setorder(
  target_modules_expression, endpoint, transition_order, regulator,
  topology_class, -supports, target
)
cross_endpoint_targets <- merge(
  target_modules[endpoint == "TF-gene", .(
    transition_order, transition, regulator, target, topology_class,
    tf_gene_supports = supports
  )],
  target_modules[endpoint == "factorized-chain", .(
    transition_order, transition, regulator, target, topology_class,
    chain_supports = supports, chain_regions = regions
  )],
  by = c(
    "transition_order", "transition", "regulator", "target", "topology_class"
  ), sort = FALSE
)
cross_endpoint_targets[, combined_supports := tf_gene_supports + chain_supports]
cross_endpoint_targets <- merge(
  cross_endpoint_targets,
  expression_transitions,
  by.x = c("transition_order", "transition", "target"),
  by.y = c("transition_order", "transition", "feature"),
  all.x = TRUE, sort = FALSE
)
setorder(
  cross_endpoint_targets, transition_order, regulator, topology_class,
  -combined_supports, target
)

pairwise_rows <- list()
pairwise_index <- 0L
for (endpoint in names(endpoint_spec)) {
  keys <- c(endpoint_spec[[endpoint]]$keys, "direction")
  folded <- fold_tables[[endpoint]]
  for (transition_order in 1:4) {
    donors <- intersect(
      unique(folded[state_order == transition_order, heldout_donor]),
      unique(folded[state_order == transition_order + 1L, heldout_donor])
    )
    for (donor in donors) {
      previous <- unique(folded[
        state_order == transition_order & heldout_donor == donor, ..keys
      ])
      next_state <- unique(folded[
        state_order == transition_order + 1L & heldout_donor == donor, ..keys
      ])
      shared <- nrow(merge(previous, next_state, by = keys))
      union <- nrow(previous) + nrow(next_state) - shared
      pairwise_index <- pairwise_index + 1L
      pairwise_rows[[pairwise_index]] <- data.table(
        endpoint = endpoint, transition_order = transition_order,
        transition = paste0(states$state[[transition_order]], " -> ", states$state[[transition_order + 1L]]),
        heldout_donor = donor, previous_edges = nrow(previous), next_edges = nrow(next_state),
        shared_signed_edges = shared, union_signed_edges = union,
        signed_jaccard = fifelse(union > 0, shared / union, NA_real_),
        turnover_fraction = fifelse(union > 0, 1 - shared / union, NA_real_)
      )
    }
  }
}
pairwise <- rbindlist(pairwise_rows)
pairwise_summary <- pairwise[, .(
  donor_folds = .N,
  median_signed_jaccard = median(signed_jaccard),
  q25_signed_jaccard = quantile(signed_jaccard, 0.25),
  q75_signed_jaccard = quantile(signed_jaccard, 0.75),
  median_turnover_fraction = median(turnover_fraction)
), by = .(endpoint, transition_order, transition)]

go_modules <- unique(target_modules[
  endpoint == "TF-gene" & topology_class %chin% c("gained", "lost"),
  .(transition_order, transition, regulator, topology_class, target)
])
universe <- sort(unique(stable_tables[["TF-gene"]]$target), method = "radix")
go_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = universe, columns = c("GOALL", "ONTOLOGYALL"), keytype = "SYMBOL"
))
go_map <- unique(as.data.table(go_map)[
  ONTOLOGYALL == "BP" & !is.na(GOALL), .(target = SYMBOL, go_id = GOALL)
])
go_sizes <- go_map[, .(universe_genes = uniqueN(target)), by = go_id][
  universe_genes >= 10L & universe_genes <= 500L
]
go_map <- go_map[go_id %in% go_sizes$go_id]
setkey(go_map, target)
module_sets <- go_modules[, .(
  genes = list(sort(unique(target), method = "radix"))
), by = .(transition_order, transition, regulator, topology_class)]
go_rows <- list()
for (module_index in seq_len(nrow(module_sets))) {
  module <- module_sets[module_index]
  genes <- module$genes[[1L]]
  annotated <- intersect(genes, go_map$target)
  if (length(annotated) < 5L) next
  hits <- go_map[J(annotated), nomatch = 0L][
    , .(overlap_genes = uniqueN(target)),
    by = go_id
  ]
  hits <- merge(hits, go_sizes, by = "go_id")
  hits[, `:=`(
    transition_order = module$transition_order, transition = module$transition,
    regulator = module$regulator, topology_class = module$topology_class,
    module_genes = length(genes), annotated_module_genes = length(annotated),
    annotated_universe_genes = uniqueN(go_map$target),
    p_value = phyper(
      overlap_genes - 1L, universe_genes, uniqueN(go_map$target) - universe_genes,
      length(annotated),
      lower.tail = FALSE
    )
  )]
  hits[, p_adjust := p.adjust(p_value, method = "BH")]
  go_rows[[length(go_rows) + 1L]] <- hits
}
go_enrichment <- if (length(go_rows)) rbindlist(go_rows) else data.table(go_id = character())
if (nrow(go_enrichment)) {
  terms <- as.data.table(AnnotationDbi::select(
    GO.db,
    keys = unique(go_enrichment$go_id), columns = "TERM", keytype = "GOID"
  ))
  setnames(terms, c("GOID", "TERM"), c("go_id", "term"))
  go_enrichment <- merge(go_enrichment, terms, by = "go_id", all.x = TRUE)
  setorder(go_enrichment, transition_order, regulator, topology_class, p_adjust, -overlap_genes)
  go_enrichment[, within_module_rank := seq_len(.N), by = .(transition_order, regulator, topology_class)]
}

fwrite(support_all, file.path(stage, "all_edge_signed_support.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(stable_all, file.path(stage, "stable_consensus_edges.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(threshold_summary, file.path(stage, "stability_threshold_sensitivity.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(transition_edges, file.path(stage, "stable_transition_edges.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(transition_summary, file.path(stage, "stable_transition_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(regulator_wide, file.path(stage, "stable_regulator_turnover.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(tf_concordance, file.path(stage, "cross_endpoint_TF_concordance.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(activity, file.path(stage, "stable_TF_activity_waves.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(activity_expression, file.path(stage, "stable_TF_activity_and_expression_waves.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(tf_candidates, file.path(stage, "candidate_TF_rewiring_with_expression.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(target_modules, file.path(stage, "stable_target_modules.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(target_modules_expression, file.path(stage, "stable_target_modules_with_expression.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(cross_endpoint_targets, file.path(stage, "cross_endpoint_target_modules.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(pairwise, file.path(stage, "donor_deletion_transition_sensitivity.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(pairwise_summary, file.path(stage, "donor_deletion_transition_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(go_enrichment, file.path(stage, "stable_rewiring_GO_BP.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")

if (!file.rename(stage, output)) stop("Unable to promote donor-stable rewiring outputs.")
cat("gse194122_donor_stable_rewiring\tPASS\t", output, "\n", sep = "")
