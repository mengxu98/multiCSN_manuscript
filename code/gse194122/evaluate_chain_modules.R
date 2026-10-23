suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)

evaluate_chain_modules <- function(root = ".", consensus_root = "results/gse194122/cell_budget", go_root = "results/gse194122/go_annotation", out = "results/gse194122/programs") {
root <- normalizePath(root, mustWork = TRUE)
consensus_root <- normalizePath(consensus_root, mustWork = TRUE)
go_root <- normalizePath(go_root, mustWork = TRUE)
dir.create(out, recursive = TRUE)
ref <- fread(file.path(root, "results/gse194122/donor_stability/stable_consensus_edges.tsv.gz"))
catalog <- fread(file.path(root, "results/gse194122/programs/module_catalog.tsv.gz"))
profiles <- fread(file.path(consensus_root, "module_state_profiles.tsv.gz"))
reference_profiles <- fread(file.path(root, "results/gse194122/programs/module_state_summary.tsv.gz"))
seeds <- 3101:3110
budgets <- setNames(lapply(seeds, function(s) fread(file.path(consensus_root, paste0("consensus_seed", s, ".tsv.gz")))), seeds)
keys <- c("regulator", "region", "target", "direction", "state_order")
counts <- rbindlist(lapply(names(budgets), function(s) {
  z <- budgets[[s]][endpoint == "factorized-chain", ..keys]
  z[, seed := as.integer(s)]
  z
}))[, .(budget_draws = uniqueN(seed)), by = keys]
chains <- merge(ref[endpoint == "factorized-chain"], counts, by = keys, all.x = TRUE)
chains[is.na(budget_draws), budget_draws := 0L]
chains[, module_id := paste0(regulator, " (", fifelse(direction > 0L, "+", "-"), ")")]
core <- chains[budget_draws >= 8L & regulator != target]
reference_peaks <- reference_profiles[, .(reference_peak = if (max(stable_chains) > 0L && sum(stable_chains == max(stable_chains)) == 1L) state_order[which.max(stable_chains)] else NA_integer_,
  RNA_peak = if (sum(median_gene_activity == max(median_gene_activity)) == 1L) state_order[which.max(median_gene_activity)] else NA_integer_,
  reference_rho = if (sd(stable_chains) > 0L && sd(median_gene_activity) > 0L) cor(stable_chains, median_gene_activity, method = "spearman") else NA_real_), by = module_id]
draw_peaks <- profiles[, .(peak = if (max(stable_chains) > 0L && sum(stable_chains == max(stable_chains)) == 1L) state_order[which.max(stable_chains)] else NA_integer_,
  rho = if (sd(support) > 0L && sd(median_gene_activity) > 0L) cor(support, median_gene_activity, method = "spearman") else NA_real_), by = .(module_id, seed)]
draw_peaks <- merge(draw_peaks, reference_peaks, by = "module_id")
stability <- draw_peaks[, .(draws = .N, finite_draws = sum(is.finite(rho)), peak_preserved = sum(!is.na(peak) & !is.na(reference_peak) & peak == reference_peak),
  median_rho = if (any(is.finite(rho))) median(rho, na.rm = TRUE) else NA_real_,
  RNA_support_positive_draws = sum(rho > 0L, na.rm = TRUE), RNA_support_negative_draws = sum(rho < 0L, na.rm = TRUE)), by = module_id]
core_counts <- core[, .(recurrent_chains = .N, recurrent_targets = uniqueN(target), recurrent_regions = uniqueN(region), recurrent_states = uniqueN(state_order)), by = module_id]
summary <- merge(catalog, stability, by = "module_id", all.x = TRUE)
summary <- merge(summary, reference_peaks, by = "module_id", all.x = TRUE)
summary <- merge(summary, core_counts, by = "module_id", all.x = TRUE)
for (nm in c("recurrent_chains", "recurrent_targets", "recurrent_regions", "recurrent_states")) set(summary, which(is.na(summary[[nm]])), nm, 0L)
summary[, primary_gate := recurrent_targets >= 10L & peak_preserved >= 8L]
setorderv(summary, c("primary_gate", "recurrent_targets", "peak_preserved", "module_id"), c(-1L, -1L, -1L, 1L))
event_keys <- c("regulator", "region", "target", "direction", "transition_order", "event")
edge_keys <- c("regulator", "region", "target", "direction")
events <- function(x) {
  rbindlist(lapply(1:4, function(t) {
    previous <- unique(x[endpoint == "factorized-chain" & state_order == t, ..edge_keys])
    following <- unique(x[endpoint == "factorized-chain" & state_order == t + 1L, ..edge_keys])
    lost <- previous[!following, on = edge_keys]
    gained <- following[!previous, on = edge_keys]
    lost[, `:=`(transition_order = t, event = "lost")]
    gained[, `:=`(transition_order = t, event = "gained")]
    rbind(lost, gained)
  }))
}
reference_events <- events(ref)
budget_events <- rbindlist(lapply(names(budgets), function(s) {z <- events(budgets[[s]]); z[, seed := as.integer(s)]; z}))
event_counts <- budget_events[, .(budget_draws = uniqueN(seed)), by = event_keys]
reference_events <- merge(reference_events, event_counts, by = event_keys, all.x = TRUE)
reference_events[is.na(budget_draws), budget_draws := 0L]
robust_events <- reference_events[budget_draws >= 8L]
eligible_draw_events <- merge(budget_events, robust_events[, ..event_keys], by = event_keys)
switches <- merge(eligible_draw_events[event == "lost"], eligible_draw_events[event == "gained"],
  by = c("regulator", "target", "direction", "transition_order", "seed"), allow.cartesian = TRUE, suffixes = c("_lost", "_gained"))
switches <- switches[, .(joint_draws = uniqueN(seed)), by = .(regulator, target, direction, transition_order, region_lost, region_gained)]
switches[, primary_switch_gate := joint_draws >= 8L]
setorderv(switches, c("primary_switch_gate", "joint_draws", "regulator", "target"), c(-1L, -1L, 1L, 1L))
go <- fread(file.path(go_root, "TERM2GENE.tsv.gz"))
go_names <- fread(file.path(go_root, "TERM2NAME.tsv"))
go_rows <- list()
for (mod in core[, .(genes = uniqueN(target)), by = module_id][genes >= 10L, module_id]) {
  genes <- unique(core[module_id == mod, target])
  backgrounds <- list(full_stable_targets = unique(ref[endpoint == "TF-gene", target]),
    module_reference_targets = unique(chains[module_id == mod, target]))
  for (bg in names(backgrounds)) {
    mapping <- go[target %chin% backgrounds[[bg]]]
    sizes <- mapping[, .(term_genes = uniqueN(target)), by = go_id][term_genes >= 10L & term_genes <= 500L]
    mapping <- mapping[go_id %chin% sizes$go_id]
    universe <- unique(mapping$target)
    members <- intersect(genes, universe)
    hits <- mapping[target %chin% members, .(overlap = uniqueN(target)), by = go_id]
    tests <- merge(sizes, hits, by = "go_id", all.x = TRUE)
    tests[is.na(overlap), overlap := 0L]
    tests[, `:=`(module_id = mod, background = bg, core_genes = length(genes), annotated_core = length(members), annotated_background = length(universe))]
    tests[, p_value := if (length(members)) phyper(overlap - 1L, term_genes, length(universe) - term_genes, length(members), lower.tail = FALSE) else 1]
    go_rows[[length(go_rows) + 1L]] <- tests
  }
}
go_tests <- rbindlist(go_rows)
if (nrow(go_tests)) {
  go_tests[, p_adjust_global := p.adjust(p_value, "BH"), by = background]
  go_tests <- merge(go_tests, go_names, by = "go_id", all.x = TRUE)
  setorder(go_tests, module_id, background, p_adjust_global, go_id)
}
fwrite(summary, file.path(out, "all_582_module_ranking.tsv"), sep = "\t")
fwrite(chains, file.path(out, "all_reference_chain_recurrence.tsv.gz"), sep = "\t")
fwrite(core, file.path(out, "all_recurrent_chains.tsv.gz"), sep = "\t")
fwrite(draw_peaks, file.path(out, "all_module_draw_peak_tests.tsv.gz"), sep = "\t")
fwrite(reference_events, file.path(out, "all_reference_path_events.tsv.gz"), sep = "\t")
fwrite(switches, file.path(out, "all_regional_substitutions.tsv.gz"), sep = "\t")
fwrite(go_tests, file.path(out, "all_GO_tests.tsv.gz"), sep = "\t")
if (nrow(go_tests)) fwrite(go_tests[p_adjust_global < 0.05], file.path(out, "GO_supported.tsv"), sep = "\t")

}

if (sys.nframe() == 0L) evaluate_chain_modules()
