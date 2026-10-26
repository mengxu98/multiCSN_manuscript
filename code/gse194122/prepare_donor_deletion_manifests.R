#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages(library(data.table))


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
if (!all(required %in% names(trajectory)) || !nrow(trajectory) ||
  anyDuplicated(trajectory$cell_id) || anyNA(trajectory[, ..required])) {
  stop("Invalid trajectory cell identifiers or metadata.")
}
states <- data.table(
  state_order = 1:5,
  state = c("HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast"),
  state_slug = c("01_HSC", "02_MK_E_prog", "03_Proerythroblast", "04_Erythroblast", "05_Normoblast")
)
if (!setequal(unique(trajectory$state), states$state)) stop("Unexpected state labels.")

tasks <- list()
task_index <- 0L
for (state_index in seq_len(nrow(states))) {
  state_name <- states$state[[state_index]]
  state_cells <- trajectory[state == state_name]
  donors <- sort(unique(state_cells$donor_id), method = "radix")
  for (donor in donors) {
    task_index <- task_index + 1L
    donor_slug <- gsub("[^A-Za-z0-9]+", "_", as.character(donor))
    task_slug <- paste0(states$state_slug[[state_index]], "__without_", donor_slug)
    retained <- state_cells[donor_id != donor, .(cell_id, state)]
    held_out <- state_cells[donor_id == donor, .N]
    if (nrow(retained) < 3L || held_out < 1L) stop("Invalid donor-deletion fold: ", task_slug)
    manifest_relative <- file.path("cells", paste0(task_slug, ".tsv.gz"))
    fwrite(
      retained[order(cell_id)], file.path(stage, manifest_relative),
      sep = "\t", quote = FALSE
    )
    tasks[[task_index]] <- data.table(
      task_order = task_index,
      state_order = states$state_order[[state_index]],
      state = state_name,
      state_slug = states$state_slug[[state_index]],
      heldout_donor = as.character(donor),
      task_slug = task_slug,
      full_cells = nrow(state_cells),
      heldout_cells = held_out,
      expected_cells = nrow(retained),
      cell_manifest = manifest_relative
    )
  }
}
task_table <- rbindlist(tasks)
fwrite(task_table, file.path(stage, "donor_deletion_plan.tsv"), sep = "\t", quote = FALSE)
publish_result_files(stage, output)
