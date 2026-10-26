#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages(library(data.table))

# ---- Parameters ----
settings <- list(
  "results/gse194122/networks",
  "results/gse194122/trajectory"
)
# ---------------------------------------------------------------
result_root <- normalizePath(settings[[1]], mustWork = TRUE)
output_dir <- settings[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

state_table <- data.table(
  state_order = 1:5,
  state = c("HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast"),
  slug = c("01_HSC", "02_MK_E_prog", "03_Proerythroblast", "04_Erythroblast", "05_Normoblast")
)
source("code/functions/network_contract.R")
invisible(lapply(file.path(result_root, state_table$slug), validate_network_result))
missing_results <- state_table[!dir.exists(file.path(result_root, slug)), slug]
if (length(missing_results)) stop("Missing state results: ", paste(missing_results, collapse = ", "))

read_endpoint <- function(slug, endpoint) {
  filename <- switch(endpoint,
    `TF-gene` = "tf_gene.tsv",
    `TF-region` = "tf_region.tsv",
    `region-gene` = "region_gene.tsv",
    `factorized-chain` = "triplets.tsv",
    `mediated-TF-gene` = "mediated_tf_gene.tsv"
  )
  path <- file.path(result_root, slug, filename)
  fread(path)
}

endpoint_rows <- vector("list", nrow(state_table))
tf_gene <- vector("list", nrow(state_table))
chains <- vector("list", nrow(state_table))

validate_ordinal <- function(direction, evidence, weight, label) {
  direction <- sign(as.numeric(direction))
  evidence <- as.numeric(evidence)
  weight <- as.numeric(weight)
  if (any(!is.finite(direction) | direction == 0) || any(!is.finite(evidence) | evidence < 0) ||
    any(!is.finite(weight))) {
    stop("Invalid weight inputs in ", label)
  }
  if (!length(weight)) {
    return(invisible(TRUE))
  }
  expected <- multiCSN::network_ordinal_weight(direction, evidence)
  if (max(abs(expected - weight)) > 1e-12) stop("Ordinal weight mismatch in ", label)
  invisible(TRUE)
}

for (i in seq_len(nrow(state_table))) {
  slug <- state_table$slug[[i]]
  endpoint_summary <- fread(file.path(result_root, slug, "endpoint_summary.tsv"))
  endpoint_rows[[i]] <- endpoint_summary[, `:=`(
    state = state_table$state[[i]], state_order = state_table$state_order[[i]]
  )]
  direct_full <- read_endpoint(slug, "TF-gene")
  tf_region_full <- read_endpoint(slug, "TF-region")
  region_gene_full <- read_endpoint(slug, "region-gene")
  chain_full <- read_endpoint(slug, "factorized-chain")
  mediated_full <- read_endpoint(slug, "mediated-TF-gene")
  validate_ordinal(direct_full$standardized_beta, direct_full$deletion_delta_bic, direct_full$weight, paste(slug, "TF-gene"))
  validate_ordinal(tf_region_full$standardized_beta, tf_region_full$deletion_delta_bic, tf_region_full$weight, paste(slug, "TF-region"))
  validate_ordinal(region_gene_full$standardized_beta, region_gene_full$deletion_delta_bic, region_gene_full$weight, paste(slug, "region-gene"))
  validate_ordinal(sign(chain_full$weight), chain_full$chain_delta_bic, chain_full$weight, paste(slug, "factorized-chain"))
  validate_ordinal(sign(mediated_full$weight), mediated_full$chain_delta_bic, mediated_full$weight, paste(slug, "mediated-TF-gene"))
  tf_gene[[i]] <- direct_full[, .(
    regulator, target,
    direction = sign(standardized_beta)
  )]
  chains[[i]] <- chain_full[, .(
    regulator, region, target,
    direction = sign(weight)
  )]
  if (anyDuplicated(tf_gene[[i]], by = c("regulator", "target"))) {
    stop("Duplicate TF-gene keys in ", slug)
  }
  if (anyDuplicated(chains[[i]], by = c("regulator", "region", "target"))) {
    stop("Duplicate chain keys in ", slug)
  }
  if (anyDuplicated(mediated_full, by = c("regulator", "target"))) {
    stop("Duplicate mediated TF-gene keys in ", slug)
  }
}
endpoint_summary <- rbindlist(endpoint_rows)
setcolorder(endpoint_summary, c("state_order", "state", "endpoint", "selected"))

classify_transition <- function(previous, next_state, keys, endpoint, transition_order) {
  merged <- multiCSN::classify_edge_transitions(previous, next_state, keys)
  merged[, `:=`(
    endpoint = endpoint,
    transition_order = transition_order,
    previous_state = state_table$state[[transition_order]],
    next_state = state_table$state[[transition_order + 1L]],
    transition = paste0(state_table$state[[transition_order]], " -> ", state_table$state[[transition_order + 1L]])
  )]
  merged[]
}

transition_edges <- list()
index <- 0L
for (transition_order in 1:4) {
  index <- index + 1L
  transition_edges[[index]] <- classify_transition(
    tf_gene[[transition_order]], tf_gene[[transition_order + 1L]],
    c("regulator", "target"), "TF-gene", transition_order
  )
  index <- index + 1L
  transition_edges[[index]] <- classify_transition(
    chains[[transition_order]], chains[[transition_order + 1L]],
    c("regulator", "region", "target"), "factorized-chain", transition_order
  )
}
transition_edges <- rbindlist(transition_edges, fill = TRUE)
transition_summary <- transition_edges[, .(edges = .N), by = .(
  endpoint, transition_order, previous_state, next_state, transition, topology_class
)]
transition_grid <- CJ(
  endpoint = c("TF-gene", "factorized-chain"),
  transition_order = 1:4,
  topology_class = c("gained", "lost", "shared_sign_stable", "shared_sign_flip"),
  unique = TRUE
)
transition_grid[, `:=`(
  previous_state = state_table$state[transition_order],
  next_state = state_table$state[transition_order + 1L]
)]
transition_grid[, transition := paste0(previous_state, " -> ", next_state)]
transition_summary <- merge(
  transition_grid, transition_summary,
  by = c("endpoint", "transition_order", "previous_state", "next_state", "transition", "topology_class"),
  all.x = TRUE, sort = FALSE
)
transition_summary[is.na(edges), edges := 0L]
transition_summary[, union_edges := sum(edges), by = .(endpoint, transition_order)]
transition_summary[, fraction := edges / union_edges]
setorder(transition_summary, endpoint, transition_order, topology_class)

regulator_summary <- transition_edges[, .(edges = .N), by = .(
  endpoint, transition_order, previous_state, next_state, transition, regulator, topology_class
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
  net_edge_change = gained - lost
)]
setorder(regulator_wide, endpoint, transition_order, -turnover, regulator)

display_regulators <- regulator_wide[turnover > 0,
  {
    threshold <- sort(unique(turnover), decreasing = TRUE)[min(8L, uniqueN(turnover))]
    .SD[turnover >= threshold]
  },
  by = .(endpoint, transition_order)
]

fwrite(endpoint_summary, file.path(output_dir, "GSE5c_state_endpoint_edge_counts.tsv"), sep = "\t")
fwrite(transition_summary, file.path(output_dir, "GSE5d_adjacent_transition_topology.tsv"), sep = "\t")
fwrite(regulator_wide, file.path(output_dir, "GSE5e_regulator_transition_counts_full.tsv.gz"), sep = "\t")
fwrite(display_regulators, file.path(output_dir, "GSE5e_regulator_transition_counts_display.tsv"), sep = "\t")
