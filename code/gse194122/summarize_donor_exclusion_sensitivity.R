suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)

summarize_donor_exclusion_sensitivity <- function(run = "results/gse194122/cell_budget/fits", reference = ".", out = "results/gse194122/cell_budget") {
run <- normalizePath(run, mustWork = TRUE)
reference <- normalizePath(reference, mustWork = TRUE)
dir.create(out, recursive = TRUE)
plan <- fread(file.path(run, "task_plan.tsv"))
stable_ref <- fread(file.path(reference, "results/gse194122/donor_stability/stable_consensus_edges.tsv.gz"))
catalog <- fread(file.path(reference, "results/gse194122/programs/module_catalog.tsv.gz"))
rna <- fread(file.path(reference, "results/gse194122/programs/module_state_summary.tsv.gz"))[, .(module_id, state_order, median_gene_activity)]
members <- unique(stable_ref[endpoint == "factorized-chain", .(regulator, direction, target)])
members[, module_id := paste0(regulator, " (", fifelse(direction > 0, "+", "-"), ")")]
members <- members[module_id %chin% catalog$module_id]
gene_universe <- sort(unique(stable_ref[endpoint == "TF-gene", target]))
go <- fread("results/gse194122/go_annotation/TERM2GENE.tsv.gz")[target %chin% gene_universe]
go_names <- fread("results/gse194122/go_annotation/TERM2NAME.tsv")
sizes <- go[, .(universe_genes = uniqueN(target)), by = go_id][universe_genes >= 10 & universe_genes <= 500]
go <- go[go_id %chin% sizes$go_id]
annotated_universe <- unique(go$target)
read_endpoint <- function(path, endpoint) {
  file <- if (endpoint == "TF-gene") "tf_gene.tsv.gz" else "triplets.tsv.gz"
  x <- fread(file.path(path, file))
  if (endpoint == "TF-gene") return(x[, .(regulator, target, direction = sign(standardized_beta))])
  x[, .(regulator, region, target, direction = chain_direction)]
}
enrich <- function(genes, seed, tf, transition, type) {
  genes <- intersect(genes, annotated_universe)
  hits <- go[target %chin% genes, .(overlap_genes = uniqueN(target)), by = go_id]
  x <- merge(sizes, hits, by = "go_id", all.x = TRUE)
  x[is.na(overlap_genes), overlap_genes := 0L]
  x[, `:=`(seed = seed, regulator = tf, transition_order = transition, topology_class = type, annotated_module_genes = length(genes), annotated_universe_genes = length(annotated_universe))]
  x[, p_value := if (length(genes)) phyper(overlap_genes - 1L, universe_genes, length(annotated_universe) - universe_genes, length(genes), lower.tail = FALSE) else 1]
  x[, p_adjust := p.adjust(p_value, "BH")]
  merge(x, go_names, by = "go_id", all.x = TRUE)
}
extract_set <- function(previous, next_state, tf, type) {
  a <- unique(previous[regulator == tf, target]); b <- unique(next_state[regulator == tf, target])
  if (type == "lost") setdiff(a, b) else setdiff(b, a)
}
cases <- data.table(regulator = c("HMGA1", "BNC2", "BNC2"), transition_order = c(4L, 3L, 4L), topology_class = c("lost", "gained", "lost"))
base_sets <- vector("list", nrow(cases)); go_rows <- list(); set_rows <- list(); coupling <- list(); profiles <- list()
for (i in seq_len(nrow(cases))) {
  item <- cases[i]
  base_sets[[i]] <- extract_set(stable_ref[endpoint == "TF-gene" & state_order == item$transition_order], stable_ref[endpoint == "TF-gene" & state_order == item$transition_order + 1L], item$regulator, item$topology_class)
  go_rows[[length(go_rows) + 1L]] <- enrich(base_sets[[i]], 0L, item$regulator, item$transition_order, item$topology_class)
}
per_seed <- function(seed_id) {
  coupling <- list(); profiles <- list(); set_rows <- list(); go_rows <- list()

  consensus <- list(); draws <- plan[seed == seed_id]
  for (state_id in 1:5) for (endpoint in c("TF-gene", "factorized-chain")) {
    tasks <- draws[state_order == state_id]
    full <- tasks[heldout_donor == "full"]
    folds <- tasks[heldout_donor != "full"]
    keys <- if (endpoint == "TF-gene") c("regulator", "target", "direction") else c("regulator", "region", "target", "direction")
    base <- read_endpoint(file.path(run, full$task_slug), endpoint)
    counted <- rbindlist(lapply(folds$task_slug, function(slug) {
      x <- read_endpoint(file.path(run, slug), endpoint)
      x[, fold := slug]
      x
    }))[, .(signed_folds = uniqueN(fold)), by = keys]
    supported <- merge(base, counted, by = keys)
    supported <- supported[signed_folds / nrow(folds) >= 0.8]
    supported[, `:=`(endpoint = endpoint, state_order = state_id, seed = seed_id, total_folds = nrow(folds))]
    consensus[[paste(state_id, endpoint)]] <- supported
  }
  all <- rbindlist(consensus, fill = TRUE)
  fwrite(all, file.path(out, paste0("consensus_seed", seed_id, ".tsv.gz")), sep = "\t")
  support <- merge(all[endpoint == "factorized-chain"], members, by = c("regulator", "direction", "target"), allow.cartesian = TRUE)
  counts <- support[, .(stable_chains = .N), by = .(module_id, state_order)]
  grid <- merge(CJ(module_id = catalog$module_id, state_order = 1:5), counts, by = c("module_id", "state_order"), all.x = TRUE)
  grid[is.na(stable_chains), stable_chains := 0L]
  grid <- merge(grid, catalog[, .(module_id, target_genes)], by = "module_id")
  grid[, support := stable_chains / target_genes]
  grid <- merge(grid, rna, by = c("module_id", "state_order"))
  grid[, seed := seed_id]
  profiles[[length(profiles) + 1L]] <- grid
  coupling[[length(coupling) + 1L]] <- grid[, .(rho = if (sd(support) > 0 && sd(median_gene_activity) > 0) cor(support, median_gene_activity, method = "spearman") else NA_real_, peak_state = if (sum(support == max(support)) == 1L) state_order[which.max(support)] else NA_integer_, total_support = sum(stable_chains)), by = .(seed, module_id)]
  for (i in seq_len(nrow(cases))) {
    item <- cases[i]
    genes <- extract_set(all[endpoint == "TF-gene" & state_order == item$transition_order], all[endpoint == "TF-gene" & state_order == item$transition_order + 1L], item$regulator, item$topology_class)
    base <- base_sets[[i]]; shared <- length(intersect(genes, base)); union <- length(union(genes, base))
    set_rows[[length(set_rows) + 1L]] <- data.table(seed = seed_id, regulator = item$regulator, transition_order = item$transition_order, topology_class = item$topology_class, baseline_genes = length(base), budget_genes = length(genes), shared_genes = shared, baseline_recall = if (length(base)) shared / length(base) else NA_real_, jaccard = if (union) shared / union else NA_real_, outside_fixed_go_universe = sum(!genes %chin% gene_universe))
    go_rows[[length(go_rows) + 1L]] <- enrich(genes, seed_id, item$regulator, item$transition_order, item$topology_class)
  }

  list(coupling = coupling, profiles = profiles, set_rows = set_rows, go_rows = go_rows)
}
seed_results <- parallel::mclapply(unique(plan$seed), per_seed, mc.cores = 10L, mc.preschedule = FALSE, mc.set.seed = FALSE)
coupling <- unlist(lapply(seed_results, `[[`, "coupling"), recursive = FALSE)
profiles <- unlist(lapply(seed_results, `[[`, "profiles"), recursive = FALSE)
set_rows <- unlist(lapply(seed_results, `[[`, "set_rows"), recursive = FALSE)
go_rows <- c(go_rows, unlist(lapply(seed_results, `[[`, "go_rows"), recursive = FALSE))

c <- rbindlist(coupling); g <- rbindlist(go_rows); s <- rbindlist(set_rows)
fwrite(c, file.path(out, "module_coupling.tsv"), sep = "\t")
fwrite(rbindlist(profiles), file.path(out, "module_state_profiles.tsv.gz"), sep = "\t")
fwrite(c[, .(modules = .N, finite_correlations = sum(is.finite(rho)), negative = sum(rho < 0, na.rm = TRUE), negative_fraction_among_finite = mean(rho[is.finite(rho)] < 0), median_rho = median(rho, na.rm = TRUE)), by = seed], file.path(out, "global_coupling_summary.tsv"), sep = "\t")
fwrite(s, file.path(out, "focal_target_set_retention.tsv"), sep = "\t")
fwrite(g, file.path(out, "focal_GO_all_tests.tsv.gz"), sep = "\t")
fwrite(g[, .(tested_terms = .N, significant_terms = sum(p_adjust < 0.05), annotated_module_genes = annotated_module_genes[[1L]]), by = .(seed, regulator, transition_order, topology_class)], file.path(out, "focal_GO_summary.tsv"), sep = "\t")
primary <- fread("results/gse194122/go_annotation/displayed_primary_terms.tsv")
keys <- c("regulator", "transition_order", "topology_class", "go_id")
primary_grid <- rbindlist(lapply(c(0L, unique(plan$seed)), function(seed_id) {x <- copy(primary); x[, seed := seed_id]; x}))
primary_tests <- merge(primary_grid, g[, .(seed, regulator, transition_order, topology_class, go_id, p_adjust, annotated_module_genes, universe_genes)], by = c("seed", keys), all.x = TRUE)
primary_tests[, `:=`(eligible_term = is.finite(p_adjust), significant = is.finite(p_adjust) & p_adjust < 0.05)]
fwrite(primary_tests, file.path(out, "displayed_GO_tests.tsv"), sep = "\t")
primary_retention <- primary_tests[seed != 0, .(supported = all(significant), eligible = all(eligible_term)), by = .(seed, regulator, go_id, term)]
primary_retention <- primary_retention[, .(draws = .N, supported_draws = sum(supported), eligible_draws = sum(eligible), support_fraction = mean(supported)), by = .(regulator, go_id, term)]
primary_baseline <- primary_tests[seed == 0, .(reference_significant = all(significant), reference_eligible = all(eligible_term), reference_worst_q = if (all(eligible_term)) max(p_adjust) else NA_real_), by = .(regulator, go_id)]
primary_retention <- merge(primary_retention, primary_baseline, by = c("regulator", "go_id"))
fwrite(primary_retention, file.path(out, "displayed_GO_reproducibility.tsv"), sep = "\t")

}

if (sys.nframe() == 0L) summarize_donor_exclusion_sensitivity()
