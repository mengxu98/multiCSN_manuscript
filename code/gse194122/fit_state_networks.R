#!/usr/bin/env Rscript
source("code/functions/network_fit.R")

fit_canonical_states <- function() {
states <- c("HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast")
slugs <- c("01_HSC", "02_MK_E_prog", "03_Proerythroblast", "04_Erythroblast", "05_Normoblast")
expected_cells <- c(1072L, 884L, 2300L, 4916L, 1780L)
cell_path <- "results/gse194122/trajectory/erythroid_cell_ids.tsv"
if (file.exists(cell_path)) stop("Output already exists: ", cell_path)
cell_rows <- vector("list", length(states))
for (i in seq_along(states)) {
  cell_rows[[i]] <- fit_state_networks(
    dataset = "GSE194122",
    input = "data/gse194122/gse194122_candidates.rds",
    output = file.path("results", "gse194122", "networks", slugs[[i]]),
    state = states[[i]],
    expected_cells = expected_cells[[i]],
    response_chunk = 1024L,
    profile = "gse194122",
    cell_manifest = NULL
  )
}
cell_table <- rbindlist(cell_rows)
if (nrow(cell_table) != sum(expected_cells) || anyDuplicated(cell_table$cell_id)) {
  stop("Fitted state cells do not define the complete erythroid universe.")
}
setorder(cell_table, cell_id)
dir.create(dirname(cell_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(cell_table, cell_path)
}

if (sys.nframe() == 0L) fit_canonical_states()
