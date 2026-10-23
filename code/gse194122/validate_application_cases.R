suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)

validate_application_cases <- function(reference = ".", consensus_root = "results/gse194122/cell_budget", go_root = "results/gse194122/go_annotation", out = "results/gse194122/cases") {
reference <- normalizePath(reference, mustWork = TRUE)
consensus_root <- normalizePath(consensus_root, mustWork = TRUE)
go_root <- normalizePath(go_root, mustWork = TRUE)
dir.create(out, recursive = TRUE)
ref <- fread(file.path(reference, "results/gse194122/donor_stability/stable_consensus_edges.tsv.gz"))
seeds <- 3101:3110
budgets <- setNames(lapply(seeds, function(s) fread(file.path(consensus_root, paste0("consensus_seed", s, ".tsv.gz")))), seeds)
go <- unique(fread(file.path(go_root, "TERM2GENE.tsv.gz")))
go_names <- fread(file.path(go_root, "TERM2NAME.tsv"))
primary <- fread(file.path(go_root, "displayed_primary_terms.tsv"))[regulator == "HMGA1"]
losses <- function(x, tf = NULL) {
  a <- unique(x[endpoint == "TF-gene" & state_order == 4L, .(regulator, target)])
  b <- unique(x[endpoint == "TF-gene" & state_order == 5L, .(regulator, target)])
  z <- a[!b, on = .(regulator, target)]
  sort(unique(if (is.null(tf)) z$target else z[regulator == tf, target]))
}
base <- losses(ref, "HMGA1")
sets <- c(list(reference = base), setNames(lapply(budgets, losses, tf = "HMGA1"), seeds))
universes <- list(full_stable_targets = sort(unique(ref[endpoint == "TF-gene", target])),
                  all_TF_transition_losses = losses(ref),
                  HMGA1_source_targets = sort(unique(ref[endpoint == "TF-gene" & regulator == "HMGA1" & state_order == 4L, target])))
retention <- rbindlist(lapply(names(sets), function(id) {
  genes <- sets[[id]]
  shared <- length(intersect(genes, base))
  data.table(set = id, baseline_genes = length(base), set_genes = length(genes), shared_genes = shared,
             baseline_recall = shared / length(base), budget_containment = shared / length(genes),
             jaccard = shared / length(union(genes, base)))
}))
fwrite(retention, file.path(out, "HMGA1_target_retention.tsv"), sep = "\t")
gene_ledger <- data.table(target = sort(unique(unlist(sets))))
gene_ledger[, in_reference := target %chin% base]
for (s in names(budgets)) gene_ledger[, (paste0("seed", s)) := target %chin% sets[[s]]]
gene_ledger[, budget_draws := rowSums(.SD), .SDcols = paste0("seed", seeds)]
fwrite(gene_ledger, file.path(out, "HMGA1_gene_recurrence.tsv"), sep = "\t")
sets$core_8_of_10 <- gene_ledger[in_reference & budget_draws >= 8L, target]
go_results <- list()
rarefaction <- list()
for (bg in names(universes)) {
  mapping <- go[target %chin% universes[[bg]]]
  sizes <- mapping[, .(term_genes = uniqueN(target)), by = go_id][term_genes >= 10L & term_genes <= 500L]
  mapping <- mapping[go_id %chin% sizes$go_id]
  universe <- unique(mapping$target)
  for (id in names(sets)) {
    genes <- intersect(sets[[id]], universe)
    hits <- mapping[target %chin% genes, .(overlap = uniqueN(target)), by = go_id]
    x <- merge(sizes, hits, by = "go_id", all.x = TRUE)
    x[is.na(overlap), overlap := 0L]
    x[, `:=`(background = bg, set = id, input_genes = length(sets[[id]]), annotated_set_genes = length(genes),
             background_genes = length(universes[[bg]]), annotated_background_genes = length(universe),
             outside_background = sum(!sets[[id]] %chin% universes[[bg]]))]
    x[, `:=`(a = overlap, b = annotated_set_genes - overlap, c = term_genes - overlap,
             d = annotated_background_genes - term_genes - annotated_set_genes + overlap)]
    x[, odds_ratio := (a * d) / (b * c)]
    x[, p_value := if (length(genes)) phyper(overlap - 1L, term_genes, length(universe) - term_genes, length(genes), lower.tail = FALSE) else 1]
    x[, p_adjust := p.adjust(p_value, "BH")]
    go_results[[length(go_results) + 1L]] <- merge(x, go_names, by = "go_id", all.x = TRUE)
  }
  baseline <- intersect(base, universe)
  for (id in as.character(seeds)) {
    shared <- intersect(sets[[id]], baseline)
    for (term_id in primary$go_id) {
      members <- mapping[go_id == term_id, target]
      K <- length(intersect(baseline, members))
      N <- length(baseline)
      n <- length(shared)
      k <- length(intersect(shared, members))
      rarefaction[[length(rarefaction) + 1L]] <- data.table(background = bg, set = id, go_id = term_id,
        eligible = term_id %chin% sizes$go_id, baseline_annotated = N, baseline_term_genes = K,
        shared_annotated = n, shared_term_genes = k, size_only_expected = n * K / N,
        size_only_lower = qhyper(0.025, K, N - K, n), size_only_upper = qhyper(0.975, K, N - K, n),
        size_only_upper_tail = phyper(k - 1L, K, N - K, n, lower.tail = FALSE))
    }
  }
}
all_go <- rbindlist(go_results)
fwrite(all_go, file.path(out, "HMGA1_GO_all_tests.tsv.gz"), sep = "\t")
display <- merge(CJ(background = names(universes), set = names(sets), go_id = primary$go_id),
                 all_go, by = c("background", "set", "go_id"), all.x = TRUE)
display <- merge(display[, !"term"], primary[, .(go_id, term)], by = "go_id", all.x = TRUE)
display[, eligible := !is.na(p_adjust)]
fwrite(display, file.path(out, "HMGA1_displayed_GO_tests.tsv"), sep = "\t")
fwrite(rbindlist(rarefaction), file.path(out, "HMGA1_size_only_reference.tsv"), sep = "\t")
chain_keys <- c("regulator", "region", "target", "direction", "state_order")
chains <- ref[endpoint == "factorized-chain" & regulator == "GATA1" & target != "GATA1"]
counts <- rbindlist(lapply(names(budgets), function(s) {
  z <- budgets[[s]][endpoint == "factorized-chain" & regulator == "GATA1", ..chain_keys]
  z[, seed := s]
  z
}))[, .(budget_draws = uniqueN(seed)), by = chain_keys]
chains <- merge(chains, counts, by = chain_keys, all.x = TRUE)
chains[is.na(budget_draws), budget_draws := 0L]
regions <- unique(chains[, .(region)])
parts <- tstrsplit(regions$region, "-")
regions[, `:=`(chrom = parts[[1L]], start = as.integer(parts[[2L]]), end = as.integer(parts[[3L]]))]
peaks <- fread(file.path(reference, "data/gse274113/external_resources/primary_erythroid_ChIP_peaks_hg38.tsv.gz"))[regulator == "GATA1"]
setkey(peaks, chrom, start, end)
overlaps <- foverlaps(regions, peaks, by.x = c("chrom", "start", "end"), by.y = c("chrom", "start", "end"), nomatch = NULL)
chip <- overlaps[, .(ChIP_peaks = .N), by = region]
chains <- merge(chains, chip, by = "region", all.x = TRUE)
chains[is.na(ChIP_peaks), ChIP_peaks := 0L]
effects <- rbindlist(lapply(c(7, 9, 11, 14), function(day) {
  x <- fread(file.path(reference, "data/gse274113/heldout_pseudobulk_full/rna_effects", paste0("day_", day, "__GATA1.tsv.gz")))
  x[, .(target = feature, day = day, logFC, FDR, effect_definition)]
}))
evaluated <- merge(chains, effects, by = "target", all.x = TRUE, allow.cartesian = TRUE)
evaluated[, `:=`(RNA_tested = is.finite(FDR) & is.finite(logFC),
                 RNA_concordant = is.finite(FDR) & FDR <= 0.05 & is.finite(logFC) & sign(logFC) == -direction)]
support <- evaluated[, .(RNA_tested_timepoints = sum(RNA_tested),
                          RNA_concordant_timepoints = sum(RNA_concordant),
                          concordant_days = paste(sort(day[RNA_concordant]), collapse = ";")), by = chain_keys]
chains <- merge(chains, support, by = chain_keys)
chains[, strict_pass := budget_draws >= 8L & ChIP_peaks > 0L & RNA_concordant_timepoints > 0L]
stages <- chains[, .(reference_chains = .N, reference_targets = uniqueN(target),
                    reference_ChIP_chains = sum(ChIP_peaks > 0L),
                    reference_RNA_concordant_chains = sum(RNA_concordant_timepoints > 0L),
                    reference_ChIP_RNA_chains = sum(ChIP_peaks > 0L & RNA_concordant_timepoints > 0L),
                    recurrent_chains = sum(budget_draws >= 8L),
                    recurrent_ChIP_chains = sum(budget_draws >= 8L & ChIP_peaks > 0L),
                    recurrent_ChIP_tested_chains = sum(budget_draws >= 8L & ChIP_peaks > 0L & RNA_tested_timepoints > 0L),
                    strict_chains = sum(strict_pass), strict_targets = uniqueN(target[strict_pass])), by = .(state_order, state)]
setorder(stages, state_order)
fwrite(chains, file.path(out, "GATA1_all_reference_chain_tests.tsv.gz"), sep = "\t")
fwrite(evaluated, file.path(out, "GATA1_all_timepoint_tests.tsv.gz"), sep = "\t")
fwrite(stages, file.path(out, "GATA1_gate_stages.tsv"), sep = "\t")
ranked <- if (any(chains$strict_pass)) chains[strict_pass == TRUE, .(passing_state_paths = .N,
                               concordant_timepoints = max(RNA_concordant_timepoints),
                               original_support = max(signed_support)), by = target] else
  data.table(target = character(), passing_state_paths = integer(), concordant_timepoints = integer(), original_support = numeric())
setorderv(ranked, c("passing_state_paths", "concordant_timepoints", "original_support", "target"), c(-1L, -1L, -1L, 1L))
fwrite(ranked, file.path(out, "GATA1_qualifying_target_ranking.tsv"), sep = "\t")
selected <- if (nrow(ranked)) ranked$target[[1L]] else character()
fwrite(evaluated[target %chin% selected], file.path(out, "GATA1_selected_locus_all_tests.tsv"), sep = "\t")

}

if (sys.nframe() == 0L) validate_application_cases()
