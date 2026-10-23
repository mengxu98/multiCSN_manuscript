#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages(library(data.table))


allocate_proportionally <- function(counts, budget) {
  exact <- counts * budget / sum(counts)
  allocated <- floor(exact)
  remainder <- budget - sum(allocated)
  if (remainder > 0L) {
    priority <- order(-(exact - allocated), names(counts), method = "radix")
    allocated[priority[seq_len(remainder)]] <- allocated[priority[seq_len(remainder)]] + 1L
  }
  if (sum(allocated) != budget || any(allocated > counts) || any(allocated < 1L)) {
    stop("Invalid proportional allocation.", call. = FALSE)
  }
  allocated
}

# ---- Parameters ----
arguments <- list(
  "results/gse194122/trajectory/erythroid_coordinates_and_pseudotime.tsv.gz",
  "results/gse194122/sampling"
)
# ---------------------------------------------------------------
trajectory_path <- arguments[[1L]]
output <- arguments[[2L]]
trajectory_path <- normalizePath(trajectory_path, winslash = "/", mustWork = TRUE)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
stage <- tempfile(".result-files-", tmpdir = dirname(output))
if (file.exists(stage)) stop("Staging path already exists: ", stage, call. = FALSE)
dir.create(stage, recursive = FALSE)
dir.create(file.path(stage, "cells"), recursive = FALSE)

trajectory <- fread(trajectory_path)
required <- c("cell_id", "state", "donor_id")
if (!all(required %in% names(trajectory)) || nrow(trajectory) != 10952L ||
    anyDuplicated(trajectory$cell_id) || anyNA(trajectory[, ..required])) {
  stop("Current trajectory table does not define the registered cell universe.")
}
states <- data.table(
  state_order = 1:5,
  state = c("HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast"),
  state_slug = c("01_HSC", "02_MK_E_prog", "03_Proerythroblast", "04_Erythroblast", "05_Normoblast")
)
if (!setequal(unique(trajectory$state), states$state)) stop("Unexpected state labels.")

budget <- 800L
seeds <- 3101:3110
tasks <- list()
allocations <- list()
task_index <- 0L
for (seed in seeds) {
  for (state_index in seq_len(nrow(states))) {
    task_index <- task_index + 1L
    state_name <- states$state[[state_index]]
    state_cells <- trajectory[state == state_name, .(cell_id, state, donor_id)][
      order(donor_id, cell_id)
    ]
    counts <- state_cells[, .N, by = donor_id][order(donor_id)]
    donor_counts <- counts$N
    names(donor_counts) <- as.character(counts$donor_id)
    donor_sample <- allocate_proportionally(donor_counts, budget)
    set.seed(seed * 100L + state_index)
    selected <- state_cells[, {
      n_take <- unname(donor_sample[as.character(donor_id[[1L]])])
      .SD[sample.int(.N, n_take)]
    }, by = donor_id][order(cell_id)]
    if (nrow(selected) != budget || anyDuplicated(selected$cell_id)) {
      stop("Invalid selected-cell manifest.")
    }
    task_slug <- sprintf("%s__budget%d_seed%d", states$state_slug[[state_index]], budget, seed)
    manifest_relative <- file.path("cells", paste0(task_slug, ".tsv.gz"))
    fwrite(selected[, .(cell_id, state)], file.path(stage, manifest_relative), sep = "\t")
    tasks[[task_index]] <- data.table(
      task_order = task_index,
      state_order = states$state_order[[state_index]],
      state = state_name,
      state_slug = states$state_slug[[state_index]],
      seed = seed,
      task_slug = task_slug,
      full_cells = nrow(state_cells),
      expected_cells = budget,
      cell_manifest = manifest_relative
    )
    allocations[[task_index]] <- data.table(
      task_slug = task_slug,
      state = state_name,
      seed = seed,
      donor_id = counts$donor_id,
      available_cells = counts$N,
      sampled_cells = unname(donor_sample),
      available_fraction = counts$N / sum(counts$N),
      sampled_fraction = unname(donor_sample) / budget
    )
  }
}
task_table <- rbindlist(tasks)
allocation_table <- rbindlist(allocations)
sampled_totals <- allocation_table[, .(sampled = sum(sampled_cells)), by = task_slug]
if (nrow(task_table) != 50L || any(task_table$expected_cells != budget) ||
    any(sampled_totals$sampled != budget) ||
    max(abs(allocation_table$sampled_fraction - allocation_table$available_fraction)) > 1 / budget) {
  stop("Cell-budget manifest self-check failed.")
}
fwrite(task_table, file.path(stage, "cell_budget_sampling_plan.tsv"), sep = "\t")
fwrite(allocation_table, file.path(stage, "donor_allocation.tsv"), sep = "\t")
publish_result_files(stage, output)
cat("gse194122_cell_budget_manifests\tPASS\t", output, "\n", sep = "")
