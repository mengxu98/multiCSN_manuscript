#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
stable_root <- file.path(repo, "results", "gse194122", "donor_stability")
expression_root <- file.path(repo, "results", "gse194122", "trajectory")
output <- file.path(repo, "results", "gse194122", "programs")
inputs <- c(
  edges = file.path(stable_root, "stable_consensus_edges.tsv.gz"),
  expression = file.path(expression_root, "donor_state_expression.tsv.gz"),
  group_counts = file.path(expression_root, "donor_state_group_counts.tsv")
)
if (any(!file.exists(inputs))) stop("Missing input: ", paste(inputs[!file.exists(inputs)], collapse = ", "))

states <- c("HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast")
focal_regulators <- c("SPI1", "GATA2", "RUNX1", "MYB", "GATA1", "TAL1", "KLF1", "NFE2")
min_module_genes <- 10L

stage <- tempfile(pattern = ".state_regulomes_", tmpdir = dirname(output))
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)

zscore <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

partial_r2 <- function(y, state, donor) {
  dat <- data.frame(y = y, state = factor(state, levels = states), donor = factor(donor))
  dat <- dat[is.finite(dat$y), , drop = FALSE]
  if (nrow(dat) < 15L || length(unique(dat$state)) < 3L || length(unique(dat$donor)) < 3L) {
    return(c(state = NA_real_, donor = NA_real_, full = NA_real_))
  }
  full <- stats::lm(y ~ state + donor, data = dat)
  donor_only <- stats::lm(y ~ donor, data = dat)
  state_only <- stats::lm(y ~ state, data = dat)
  sse_full <- sum(stats::residuals(full)^2)
  sse_donor <- sum(stats::residuals(donor_only)^2)
  sse_state <- sum(stats::residuals(state_only)^2)
  sst <- sum((dat$y - mean(dat$y))^2)
  c(
    state = if (sse_donor > 0) (sse_donor - sse_full) / sse_donor else NA_real_,
    donor = if (sse_state > 0) (sse_state - sse_full) / sse_state else NA_real_,
    full = if (sst > 0) 1 - sse_full / sst else NA_real_
  )
}

edges <- fread(inputs[["edges"]])
edges <- edges[endpoint == "factorized-chain" & stable_consensus == TRUE & is.finite(direction)]
if (!nrow(edges) || !all(edges$state %chin% states)) stop("Stable factorized-chain input changed.")
edges[, direction_label := fifelse(direction > 0, "positive", "negative")]
edges[, module_id := paste0(regulator, " (", fifelse(direction > 0, "+", "-"), ")")]

catalog <- edges[, .(
  regulator = regulator[1L],
  direction = direction[1L],
  direction_label = direction_label[1L],
  target_genes = uniqueN(target),
  regulatory_regions = uniqueN(region),
  module_stable_chains = .N,
  represented_states = uniqueN(state),
  median_signed_support = median(signed_support, na.rm = TRUE)
), by = module_id]
catalog <- catalog[target_genes >= min_module_genes]
if (nrow(catalog) < 50L) stop("Too few donor-stable regulons after the prespecified size gate.")

membership <- unique(edges[module_id %chin% catalog$module_id, .(module_id, feature = target)])
expression <- fread(inputs[["expression"]])
groups <- fread(inputs[["group_counts"]])
groups[, state := factor(state, levels = states)]
setorder(groups, state, donor_id)
if (nrow(groups) != 49L || uniqueN(groups$donor_id) != 10L) stop("Donor-state group schema changed.")

expression <- expression[feature %chin% membership$feature]
expression[, group_key := paste(state, donor_id, sep = "\r")]
group_keys <- paste(as.character(groups$state), groups$donor_id, sep = "\r")
expr_wide <- dcast(
  expression,
  feature ~ group_key,
  value.var = "mean_normalized_expression",
  fill = 0
)
features <- expr_wide$feature
expr_mat <- as.matrix(expr_wide[, ..group_keys])
rownames(expr_mat) <- features
expr_z <- t(apply(expr_mat, 1L, zscore))

membership <- membership[feature %chin% features]
module_levels <- catalog$module_id
module_matrix <- sparseMatrix(
  i = match(membership$module_id, module_levels),
  j = match(membership$feature, features),
  x = 1,
  dims = c(length(module_levels), length(features)),
  dimnames = list(module_levels, features)
)
module_sizes_observed <- Matrix::rowSums(module_matrix)
module_score <- (module_matrix %*% expr_z) / pmax(1, module_sizes_observed)
colnames(module_score) <- group_keys

activity <- as.data.table(as.table(as.matrix(module_score)))
setnames(activity, c("module_id", "group_key", "gene_module_score"))
activity <- merge(activity, groups[, .(group_key = paste(as.character(state), donor_id, sep = "\r"), state, donor_id, cells)], by = "group_key")
activity[, state_order := match(as.character(state), states)]

regulatory <- edges[module_id %chin% catalog$module_id, .(
  stable_chains = .N,
  stable_targets = uniqueN(target),
  stable_regions = uniqueN(region),
  median_signed_support = median(signed_support, na.rm = TRUE)
), by = .(module_id, state)]
reg_grid <- CJ(module_id = module_levels, state = states, unique = TRUE)
regulatory <- merge(reg_grid, regulatory, by = c("module_id", "state"), all.x = TRUE)
for (field in c("stable_chains", "stable_targets", "stable_regions")) set(regulatory, which(is.na(regulatory[[field]])), field, 0)
regulatory[, state_order := match(state, states)]
regulatory[, regulatory_activity := stable_chains / catalog$target_genes[match(module_id, catalog$module_id)]]
regulatory[, regulatory_activity_z := zscore(regulatory_activity), by = module_id]

activity <- merge(activity, regulatory[, .(module_id, state, regulatory_activity, regulatory_activity_z, stable_chains, stable_targets, stable_regions)], by = c("module_id", "state"), all.x = TRUE)
activity <- merge(activity, catalog, by = "module_id", all.x = TRUE)
setorder(activity, module_id, state_order, donor_id)

state_summary <- activity[, .(
  donors = .N,
  cells = sum(cells),
  median_gene_activity = median(gene_module_score),
  q25_gene_activity = quantile(gene_module_score, 0.25),
  q75_gene_activity = quantile(gene_module_score, 0.75),
  regulatory_activity = regulatory_activity[1L],
  regulatory_activity_z = regulatory_activity_z[1L],
  stable_chains = stable_chains[1L],
  stable_targets = stable_targets[1L],
  stable_regions = stable_regions[1L]
), by = .(module_id, regulator, direction, direction_label, target_genes, regulatory_regions, state_order, state)]

model_rows <- activity[,
  {
    fit <- partial_r2(gene_module_score, as.character(state), donor_id)
    .(partial_state_r2 = fit[["state"]], partial_donor_r2 = fit[["donor"]], full_model_r2 = fit[["full"]])
  },
  by = .(module_id, regulator, direction, direction_label, target_genes, regulatory_regions)
]
coupling <- state_summary[, .(
  activity_coupling_rho = suppressWarnings(cor(median_gene_activity, regulatory_activity, method = "spearman")),
  peak_gene_state_order = state_order[which.max(median_gene_activity)],
  peak_regulatory_state_order = state_order[which.max(regulatory_activity)]
), by = module_id]
model_rows <- merge(model_rows, coupling, by = "module_id", all.x = TRUE)

dynamics <- copy(state_summary)
setorder(dynamics, module_id, state_order)
dynamics[, `:=`(
  delta_gene_activity = median_gene_activity - shift(median_gene_activity),
  delta_regulatory_activity = regulatory_activity_z - shift(regulatory_activity_z),
  previous_state = shift(as.character(state))
), by = module_id]
dynamics <- dynamics[state_order > 1L]
dynamics[, transition := paste(previous_state, as.character(state), sep = " -> ")]
dynamics[, direction_concordant := sign(delta_gene_activity) == sign(delta_regulatory_activity)]

eligible_hero <- merge(model_rows, catalog, by = c("module_id", "regulator", "direction", "direction_label", "target_genes", "regulatory_regions"))
eligible_hero <- eligible_hero[
  regulator %chin% focal_regulators & represented_states >= 3L &
    is.finite(partial_state_r2) & is.finite(activity_coupling_rho)
]
eligible_hero[, hero_score := frank(-partial_state_r2, ties.method = "average") + frank(-activity_coupling_rho, ties.method = "average") + 0.25 * frank(-target_genes, ties.method = "average")]
setorder(eligible_hero, hero_score, -partial_state_r2, -activity_coupling_rho, -target_genes, regulator)
heroes <- eligible_hero[, .SD[1L], by = regulator][1:min(6L, uniqueN(eligible_hero$regulator))]
if (nrow(heroes) < 4L) stop("Fewer than four prespecified focal regulons passed the hero gate.")
heroes[, hero_rank := seq_len(.N)]

hero_network_edges <- edges[module_id %chin% heroes$module_id, .(
  represented_states = uniqueN(state),
  exact_regions = uniqueN(region),
  stable_chains = .N,
  median_signed_support = median(signed_support, na.rm = TRUE)
), by = .(module_id, regulator, direction, direction_label, target)]
target_state_expression <- expression[feature %chin% hero_network_edges$target, .(
  median_expression = median(mean_normalized_expression)
), by = .(target = feature, state)]
setorder(target_state_expression, target, -median_expression, state)
target_peaks <- target_state_expression[, .(
  peak_expression_state = as.character(state[1L]),
  peak_expression = median_expression[1L]
), by = target]
hero_network_edges <- merge(hero_network_edges, target_peaks, by = "target", all.x = TRUE)
hero_network_edges[, hero_order := match(regulator, heroes$regulator)]
setorder(hero_network_edges, hero_order, -represented_states, -stable_chains, -exact_regions, target)
hero_network_edges[, hero_order := NULL]

hero_late <- heroes[regulator %chin% c("GATA1", "TAL1", "KLF1", "NFE2")]
if (!nrow(hero_late)) hero_late <- heroes
setorder(hero_late, -partial_state_r2, -activity_coupling_rho, -target_genes)
hero_module <- hero_late$module_id[[1L]]
hero_edges <- edges[module_id == hero_module]
hero_targets <- hero_edges[, .(
  represented_states = uniqueN(state),
  exact_regions = uniqueN(region),
  stable_chains = .N,
  median_signed_support = median(signed_support)
), by = .(module_id, regulator, direction, direction_label, target)]
setorder(hero_targets, -represented_states, -stable_chains, -exact_regions, -median_signed_support, target)
hero_targets <- hero_targets[1:min(12L, .N)]
hero_target_state <- hero_edges[target %chin% hero_targets$target, .(
  exact_regions = uniqueN(region),
  stable_chains = .N,
  median_signed_support = median(signed_support)
), by = .(module_id, regulator, direction, direction_label, target, state)]
hero_expr <- expression[feature %chin% hero_targets$target, .(
  median_expression = median(mean_normalized_expression),
  q25_expression = quantile(mean_normalized_expression, 0.25),
  q75_expression = quantile(mean_normalized_expression, 0.75),
  donors = .N
), by = .(target = feature, state)]
hero_target_state <- merge(CJ(target = hero_targets$target, state = states, unique = TRUE), hero_target_state, by = c("target", "state"), all.x = TRUE)
hero_target_state <- merge(hero_target_state, hero_expr, by = c("target", "state"), all.x = TRUE)
hero_target_state <- merge(hero_target_state, hero_targets, by = "target", all.x = TRUE, suffixes = c("_state", "_module"))
hero_target_state[, state_order := match(state, states)]
for (field in c("exact_regions_state", "stable_chains_state")) set(hero_target_state, which(is.na(hero_target_state[[field]])), field, 0)

fwrite(catalog, file.path(stage, "module_catalog.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(activity, file.path(stage, "module_activity_by_donor_state.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(state_summary, file.path(stage, "module_state_summary.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(model_rows, file.path(stage, "module_model_variance.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(dynamics, file.path(stage, "module_state_dynamics.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(heroes, file.path(stage, "hero_regulons.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(hero_network_edges, file.path(stage, "hero_regulon_network.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(hero_target_state, file.path(stage, "hero_regulon_targets.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")

publish_result_files(stage, output)
cat("gse194122_state_regulomes\tPASS\t", nrow(catalog), " modules; hero ", hero_module, "\n", sep = "")
