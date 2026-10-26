fit_cell_budget_networks <- function(config, root) {
root <- normalizePath(root, mustWork = FALSE)
suppressPackageStartupMessages({library(data.table); library(Matrix); library(SeuratObject); library(Signac); library(multiCSN)})
setDTthreads(1L)
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
allocate <- function(counts, budget) {
  exact <- counts * budget / sum(counts)
  n <- floor(exact)
  left <- budget - sum(n)
  if (left > 0) {
    order <- order(-(exact - n), names(counts), method = "radix")
    n[order[seq_len(left)]] <- n[order[seq_len(left)]] + 1L
  }
  n
}
dir.create(root, recursive = TRUE, showWarnings = FALSE)
lock <- file.path(root, "RUNNING")
if (!dir.create(lock, showWarnings = FALSE)) stop("A run lock already exists")
on.exit(unlink(lock, recursive = TRUE), add = TRUE)
input <- normalizePath(config$input, mustWork = TRUE)
object <- readRDS(input)
meta <- as.data.table(object[[]], keep.rownames = "cell_id")
if (config$dataset == "erythroid") {
  trajectory <- fread(config$trajectory)
  groups <- config$states
  all_tasks <- list()
  for (seed in config$seeds) for (state_index in seq_along(groups)) {
    state_name <- groups[[state_index]]
    available <- trajectory[state == state_name, .(cell_id, donor_id)][order(donor_id, cell_id)]
    folds <- c("full", sort(unique(as.character(available$donor_id))))
    for (fold in folds) {
      pool <- if (fold == "full") copy(available) else available[as.character(donor_id) != fold]
      counts <- pool[, .N, by = donor_id][order(donor_id)]
      n <- allocate(setNames(counts$N, as.character(counts$donor_id)), config$budget)
      set.seed(seed * 100L + state_index)
      selected <- pool[, .SD[sample.int(.N, n[as.character(donor_id[[1L]])])], by = donor_id][order(cell_id)]
      slug <- sprintf("seed%d_state%d_%s", seed, state_index, fold)
      all_tasks[[length(all_tasks) + 1L]] <- list(slug = slug, seed = seed, group = state_name, state_order = state_index, fold = fold, cells = selected$cell_id)
    }
  }
} else {
  class_map <- c(GluN2 = "GluN", GluN3 = "GluN", GluN4 = "GluN", GluN5 = "GluN", `nIPC/GluN1` = "GluN", SP = "GluN", IN1 = "IN", IN2 = "IN", IN3 = "IN", RG = "IPC_RG", `Cyc. Prog.` = "IPC_RG", `mGPC/OPC` = "OPC_Glia", `EC/Peric.` = "Vascular")
  meta[, broad_class := unname(class_map[author_cluster_name])]
  all_tasks <- lapply(config$classes, function(class) list(slug = class, seed = 1L, group = class, state_order = match(class, config$classes), fold = "full", cells = sort(meta[broad_class == class, cell_id])))
}
plan <- rbindlist(lapply(all_tasks, function(t) data.table(task_slug = t$slug, seed = t$seed, state = t$group, state_order = t$state_order, heldout_donor = t$fold, cells = length(t$cells))))
fwrite(plan, file.path(root, "task_plan.tsv"), sep = "\t")
saveRDS(all_tasks, file.path(root, "cell_manifests.rds"))
worker <- function(task) {
  path <- file.path(root, task$slug)
  required <- c("tf_gene.tsv.gz", "tf_region.tsv.gz", "region_gene.tsv.gz", "triplets.tsv.gz", "mediated_tf_gene.tsv.gz", "fit_summary.rds")
  if (all(file.exists(file.path(path, required)))) return(TRUE)
    dir.create(path, showWarnings = FALSE)
    fit <- multiCSN::fit_layered_network(object, cell_group = "aggregate", cells = task$cells, min_detected = 0L, min_peak_cells = 3L, upstream = 100000, downstream = 0, tf_region_scope = "domain_peaks", exclude_response_alias = FALSE, renormalize = FALSE, response_chunk = 1024L, cores = config$fit_cores, sort_regulators = TRUE, verbose = TRUE)
    outputs <- list(tf_gene = fit$tf_gene, tf_region = fit$tf_region, region_gene = fit$region_gene, triplets = fit$triplets, mediated_tf_gene = fit$mediated_tf_gene)
    for (name in names(outputs)) fwrite(outputs[[name]], file.path(path, paste0(name, ".tsv.gz")), sep = "\t")
    saveRDS(list(settings = fit$settings, summary = fit$summary, accounting = fit$accounting), file.path(path, "fit_summary.rds"))
    TRUE
}
invisible(parallel::mclapply(all_tasks, worker, mc.cores = config$parallel_jobs, mc.preschedule = FALSE, mc.set.seed = FALSE))
}

run_cortex_sensitivity <- function() {
  fit_cell_budget_networks(
    list(input = "data/gse162170/gse162170_candidates.rds", dataset = "cortex",
         classes = c("GluN", "IN", "IPC_RG", "OPC_Glia"),
         fit_cores = 1L, parallel_jobs = 2L),
    "results/gse162170/within_class/fits"
  )
}

run_erythroid_sensitivity <- function() {
  fit_cell_budget_networks(
    list(input = "data/gse194122/gse194122_candidates.rds", dataset = "erythroid",
         trajectory = "results/gse194122/trajectory/erythroid_coordinates_and_pseudotime.tsv.gz",
         states = c("HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast"),
         budget = 500L, seeds = 3101:3110, fit_cores = 1L,
         parallel_jobs = 24L),
    "results/gse194122/cell_budget/fits"
  )
}
