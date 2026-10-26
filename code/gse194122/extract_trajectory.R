#!/usr/bin/env Rscript

source("code/functions/table_io.R")

source("code/functions/hdf5_io.R")

suppressPackageStartupMessages({
  library(data.table)
  library(rhdf5)
})


orient_embedding <- function(x, observations) {
  if (identical(dim(x), c(2L, observations))) {
    return(t(x))
  }
  if (identical(dim(x), c(observations, 2L))) {
    return(x)
  }
  stop("Unexpected embedding dimensions: ", paste(dim(x), collapse = " x "))
}

# ---- Parameters ----
arguments <- list(
  "data/gse194122/processed_multiome.h5ad",
  "results/gse194122/trajectory/erythroid_cell_ids.tsv",
  "results/gse194122/trajectory"
)
# ---------------------------------------------------------------
input <- arguments[[1L]]
cell_path <- arguments[[2L]]
output <- arguments[[3L]]
input <- normalizePath(input, winslash = "/", mustWork = TRUE)
cell_path <- normalizePath(cell_path, winslash = "/", mustWork = TRUE)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
stage <- tempfile(".result-files-", tmpdir = dirname(output))
if (file.exists(stage)) stop("Staging path already exists: ", stage, call. = FALSE)
dir.create(stage, recursive = FALSE)

cell_ids <- h5read(input, "/obs/_index")
if (!length(cell_ids) || anyDuplicated(cell_ids)) {
  stop("Invalid H5AD cell identifiers.")
}
current <- fread(cell_path)
if (!identical(names(current), c("cell_id", "state")) ||
  !nrow(current) || anyDuplicated(current$cell_id)) {
  stop("Invalid erythroid cell manifest.")
}
index <- match(current$cell_id, cell_ids)
if (anyNA(index) || anyDuplicated(index)) {
  stop("Invalid cell mapping.")
}

states <- read_hdf5_category(input, "cell_type")[index]
if (anyNA(states) || any(states != current$state)) {
  stop("Cell states do not match.")
}
gex_umap <- orient_embedding(h5read(input, "/obsm/GEX_X_umap"), length(cell_ids))[index, ]
atac_umap <- orient_embedding(h5read(input, "/obsm/ATAC_umap"), length(cell_ids))[index, ]
gex_pseudotime <- h5read(input, "/obs/GEX_pseudotime_order")[index]
atac_pseudotime <- h5read(input, "/obs/ATAC_pseudotime_order")[index]

if (any(!is.finite(gex_umap)) || any(!is.finite(atac_umap)) ||
  any(!is.finite(atac_pseudotime)) ||
  any(atac_pseudotime < 0 | atac_pseudotime > 1)) {
  stop("Invalid embedding or ATAC pseudotime.")
}

trajectory <- data.table(
  cell_id = current$cell_id,
  state = current$state,
  donor_id = h5read(input, "/obs/DonorID")[index],
  site = read_hdf5_category(input, "Site")[index],
  batch = read_hdf5_category(input, "batch")[index],
  sample = read_hdf5_category(input, "Samplename")[index],
  GEX_UMAP_1 = gex_umap[, 1L],
  GEX_UMAP_2 = gex_umap[, 2L],
  ATAC_UMAP_1 = atac_umap[, 1L],
  ATAC_UMAP_2 = atac_umap[, 2L],
  GEX_pseudotime = gex_pseudotime,
  ATAC_pseudotime = atac_pseudotime
)
setorder(trajectory, cell_id)

fwrite(
  trajectory,
  file.path(stage, "erythroid_coordinates_and_pseudotime.tsv.gz"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  trajectory[, .(
    cells = .N,
    donors = uniqueN(donor_id),
    samples = uniqueN(sample),
    GEX_pseudotime_finite = sum(is.finite(GEX_pseudotime)),
    ATAC_pseudotime_finite = sum(is.finite(ATAC_pseudotime)),
    ATAC_pseudotime_median = median(ATAC_pseudotime),
    ATAC_pseudotime_q25 = quantile(ATAC_pseudotime, 0.25),
    ATAC_pseudotime_q75 = quantile(ATAC_pseudotime, 0.75)
  ), by = state],
  file.path(stage, "state_trajectory_summary.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  trajectory[, .(cells = .N), by = .(donor_id, site, batch, sample, state)],
  file.path(stage, "donor_state_cell_counts.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
publish_result_files(stage, output)
