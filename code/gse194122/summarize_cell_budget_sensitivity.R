#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)

overlap_summary <- function(reference, sampled, keys) {
  reference <- copy(reference)
  sampled <- copy(sampled)
  reference[, reference_present := TRUE]
  sampled[, sampled_present := TRUE]
  merged <- merge(reference, sampled, by = keys, all = TRUE, suffixes = c("_reference", "_sampled"))
  merged[is.na(reference_present), reference_present := FALSE]
  merged[is.na(sampled_present), sampled_present := FALSE]
  shared_unsigned <- merged[reference_present & sampled_present, .N]
  shared_signed <- merged[
    reference_present & sampled_present & direction_reference == direction_sampled,
    .N
  ]
  data.table(
    reference_edges = nrow(reference), sampled_edges = nrow(sampled),
    shared_unsigned = shared_unsigned, shared_signed = shared_signed,
    unsigned_jaccard = shared_unsigned / nrow(merged),
    signed_jaccard = shared_signed / (nrow(reference) + nrow(sampled) - shared_signed),
    sign_concordance_shared = shared_signed / shared_unsigned
  )
}

composition_summary <- function(reference, sampled, keys) {
  classes <- c("gained", "lost", "shared_sign_stable", "shared_sign_flip")
  count_classes <- function(x, prefix) {
    result <- x[, .(edges = .N), by = topology_class]
    result <- merge(data.table(topology_class = classes), result, by = "topology_class", all.x = TRUE)
    result[is.na(edges), edges := 0L]
    result[, fraction := edges / sum(edges)]
    setnames(result, c("edges", "fraction"), paste0(prefix, c("_edges", "_fraction")))
    result
  }
  fractions <- merge(
    count_classes(reference, "reference"), count_classes(sampled, "sampled"),
    by = "topology_class", sort = FALSE
  )
  joined <- merge(
    reference[, c(keys, "topology_class"), with = FALSE],
    sampled[, c(keys, "topology_class"), with = FALSE],
    by = keys, all = TRUE, suffixes = c("_reference", "_sampled")
  )
  reference_mask <- !is.na(joined$topology_class_reference)
  sampled_mask <- !is.na(joined$topology_class_sampled)
  summary <- data.table(
    composition_similarity = 1 - 0.5 * sum(abs(fractions$reference_fraction - fractions$sampled_fraction)),
    reference_key_recall = mean(sampled_mask[reference_mask]),
    exact_class_recall = mean(
      joined$topology_class_reference[reference_mask] == joined$topology_class_sampled[reference_mask],
      na.rm = TRUE
    ),
    class_agreement_shared = mean(
      joined$topology_class_reference[reference_mask & sampled_mask] ==
        joined$topology_class_sampled[reference_mask & sampled_mask]
    )
  )
  list(fractions = fractions, summary = summary)
}

regulator_turnover <- function(edges) {
  edges[, .(
    turnover = sum(topology_class %chin% c("gained", "lost", "shared_sign_flip"))
  ), by = regulator]
}

turnover_spearman <- function(reference, sampled) {
  joined <- merge(
    regulator_turnover(reference), regulator_turnover(sampled),
    by = "regulator", all = TRUE, suffixes = c("_reference", "_sampled")
  )
  joined[is.na(turnover_reference), turnover_reference := 0L]
  joined[is.na(turnover_sampled), turnover_sampled := 0L]
  if (nrow(joined) < 2L || !isTRUE(sd(joined$turnover_reference) > 0) ||
      !isTRUE(sd(joined$turnover_sampled) > 0)) return(NA_real_)
  cor(joined$turnover_reference, joined$turnover_sampled, method = "spearman")
}

# ---- Parameters ----
settings <- list(
  "runs/gse194122/cell_budget_runs/results",
  "results/gse194122/sampling/cell_budget_sampling_plan.tsv",
  "results/gse194122/networks",
  "results/gse194122/programs/hero_regulons.tsv",
  "results/gse194122/donor_stability/stable_consensus_edges.tsv.gz",
  "results/gse194122/cell_budget"
)
# ---------------------------------------------------------------
sample_root <- normalizePath(settings[[1]], mustWork = TRUE)
task_plan_path <- normalizePath(settings[[2]], mustWork = TRUE)
formal_root <- normalizePath(settings[[3]], mustWork = TRUE)
hero_path <- normalizePath(settings[[4]], mustWork = TRUE)
stable_path <- normalizePath(settings[[5]], mustWork = TRUE)
output_dir <- settings[[6]]
source("code/functions/network_contract.R")
stage <- tempfile(".result-files-", tmpdir = dirname(output_dir))
if (file.exists(stage)) stop("Staging output exists: ", stage)
dir.create(stage, recursive = TRUE, showWarnings = FALSE)

tasks <- fread(task_plan_path)
invisible(lapply(file.path(sample_root, tasks$task_slug), validate_network_result))
state_table <- unique(tasks[, .(state_order, state, state_slug)])
setorder(state_table, state_order)
endpoint_specs <- list(
  `TF-gene` = list(file = "tf_gene.tsv", keys = c("regulator", "target"), direction = "standardized_beta"),
  `factorized-chain` = list(file = "triplets.tsv", keys = c("regulator", "region", "target"), direction = "chain_direction")
)

read_endpoint <- function(slug, root, endpoint) {
  spec <- endpoint_specs[[endpoint]]
  path <- file.path(root, slug, spec$file)
  if (!file.exists(path)) stop("Missing endpoint file: ", path)
  x <- fread(path, select = c(spec$keys, spec$direction))
  setnames(x, spec$direction, "direction")
  x[, direction := sign(as.numeric(direction))]
  if (any(!is.finite(x$direction) | x$direction == 0) || anyDuplicated(x, by = spec$keys)) {
    stop("Invalid or duplicate endpoint keys in: ", path)
  }
  x
}

formal <- lapply(names(endpoint_specs), function(endpoint) {
  setNames(lapply(state_table$state_slug, read_endpoint, root = formal_root, endpoint = endpoint), state_table$state_slug)
})
names(formal) <- names(endpoint_specs)

state_overlap <- list()
transition_fractions <- list()
transition_summary <- list()
module_support <- list()
index_overlap <- index_fraction <- index_transition <- index_module <- 0L

heroes <- fread(hero_path)[, .(
  module_id, regulator, direction = as.integer(direction)
)]
stable_edges <- fread(stable_path)
hero_targets <- unique(merge(
  stable_edges[
    endpoint == "factorized-chain" & stable_consensus,
    .(regulator, direction = as.integer(direction), target)
  ],
  heroes,
  by = c("regulator", "direction"), sort = FALSE
)[, .(module_id, regulator, direction, target)])
hero_modules <- hero_targets[, .(module_targets = uniqueN(target)), by = .(module_id, regulator, direction)]

summarize_modules <- function(chains, fit_type, seed, state_row) {
  hits <- merge(
    chains, hero_targets,
    by = c("regulator", "target", "direction"), allow.cartesian = TRUE, sort = FALSE
  )
  observed <- hits[, .(
    selected_targets = uniqueN(target), selected_chains = uniqueN(paste(regulator, region, target, sep = "\r"))
  ), by = .(module_id, regulator, direction)]
  result <- merge(hero_modules, observed, by = c("module_id", "regulator", "direction"), all.x = TRUE)
  result[is.na(selected_targets), `:=`(selected_targets = 0L, selected_chains = 0L)]
  result[, `:=`(
    fit_type = fit_type, seed = seed,
    state_order = state_row$state_order, state = state_row$state,
    target_support_fraction = selected_targets / module_targets
  )]
  result[]
}

for (i in seq_len(nrow(state_table))) {
  index_module <- index_module + 1L
  module_support[[index_module]] <- summarize_modules(
    formal[["factorized-chain"]][[state_table$state_slug[[i]]]],
    "all-cell", NA_integer_, state_table[i]
  )
}

for (seed_value in sort(unique(tasks$seed))) {
  seed_tasks <- tasks[seed == seed_value][order(state_order)]
  sampled <- lapply(names(endpoint_specs), function(endpoint) {
    setNames(lapply(seed_tasks$task_slug, read_endpoint, root = sample_root, endpoint = endpoint), seed_tasks$state_slug)
  })
  names(sampled) <- names(endpoint_specs)

  for (endpoint in names(endpoint_specs)) {
    keys <- endpoint_specs[[endpoint]]$keys
    for (i in seq_len(nrow(state_table))) {
      slug <- state_table$state_slug[[i]]
      index_overlap <- index_overlap + 1L
      state_overlap[[index_overlap]] <- cbind(
        data.table(
          seed = seed_value, endpoint = endpoint,
          state_order = state_table$state_order[[i]], state = state_table$state[[i]]
        ),
        overlap_summary(formal[[endpoint]][[slug]], sampled[[endpoint]][[slug]], keys)
      )
      if (endpoint == "factorized-chain") {
        index_module <- index_module + 1L
        module_support[[index_module]] <- summarize_modules(
          sampled[[endpoint]][[slug]], "budget-800", seed_value, state_table[i]
        )
      }
    }

    for (transition_order in 1:4) {
      previous_slug <- state_table$state_slug[[transition_order]]
      next_slug <- state_table$state_slug[[transition_order + 1L]]
      reference_transition <- multiCSN::classify_edge_transitions(
        formal[[endpoint]][[previous_slug]], formal[[endpoint]][[next_slug]], keys
      )
      sampled_transition <- multiCSN::classify_edge_transitions(
        sampled[[endpoint]][[previous_slug]], sampled[[endpoint]][[next_slug]], keys
      )
      composition <- composition_summary(reference_transition, sampled_transition, keys)
      transition_label <- paste(state_table$state[[transition_order]], state_table$state[[transition_order + 1L]], sep = " -> ")
      index_fraction <- index_fraction + 1L
      transition_fractions[[index_fraction]] <- composition$fractions[, `:=`(
        seed = seed_value, endpoint = endpoint, transition_order = transition_order,
        transition = transition_label
      )]
      index_transition <- index_transition + 1L
      transition_summary[[index_transition]] <- composition$summary[, `:=`(
        seed = seed_value, endpoint = endpoint, transition_order = transition_order,
        transition = transition_label,
        regulator_turnover_spearman = turnover_spearman(reference_transition, sampled_transition)
      )]
    }
  }
}

outputs <- list(
  state_network_overlap = rbindlist(state_overlap),
  transition_topology_fractions = rbindlist(transition_fractions),
  transition_topology_summary = rbindlist(transition_summary),
  focal_module_support = rbindlist(module_support, fill = TRUE)
)

for (name in names(outputs)) {
  fwrite(outputs[[name]], file.path(stage, paste0(name, ".tsv")), sep = "\t", quote = FALSE, na = "NA")
}

publish_result_files(stage, output_dir)
cat("GSE194122_SENSITIVITY_SUMMARY_PASS\n")
