#!/usr/bin/env Rscript

source("code/functions/hdf5_io.R")

suppressPackageStartupMessages({
  library(data.table)
  library(HDF5Array)
  library(Matrix)
  library(rhdf5)
  library(Seurat)
  library(harmony)
})


# ---- Parameters ----
arguments <- list(
  "data/gse194122/processed_multiome.h5ad",
  "results/gse194122/state_transitions/erythroid_cell_ids.tsv",
  "results/gse194122/erythroid_gex_umap"
)
# ---------------------------------------------------------------
input <- arguments[[1L]]
cell_path <- arguments[[2L]]
output <- arguments[[3L]]
input <- normalizePath(input, winslash = "/", mustWork = TRUE)
cell_path <- normalizePath(cell_path, winslash = "/", mustWork = TRUE)
if (dir.exists(output)) stop("Output already exists: ", output, call. = FALSE)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
stage <- paste0(output, ".staging")
if (file.exists(stage)) stop("Staging path already exists: ", stage, call. = FALSE)
dir.create(stage, recursive = FALSE)

set.seed(20260822L)
cell_ids <- h5read(input, "/obs/_index")
if (length(cell_ids) != 69249L || anyDuplicated(cell_ids)) {
  stop("The H5AD cell universe is not the registered 69,249-cell release.")
}
current <- fread(cell_path)
if (!identical(names(current), c("cell_id", "state")) ||
  nrow(current) != 10952L || anyDuplicated(current$cell_id)) {
  stop("Current erythroid cell manifest is invalid.")
}
cell_index <- match(current$cell_id, cell_ids)
if (anyNA(cell_index) || anyDuplicated(cell_index)) {
  stop("Current erythroid cells do not map one-to-one to the H5AD.")
}

feature_names <- h5read(input, "/var/_index")
feature_codes <- h5read(input, "/var/feature_types")
feature_levels <- h5read(input, "/var/__categories/feature_types")
feature_types <- feature_levels[feature_codes + 1L]
gex_index <- which(feature_types == "GEX")
gex_names <- as.character(feature_names[gex_index])
if (length(gex_index) != 13431L || anyDuplicated(gex_names)) {
  stop("The registered 13,431-feature GEX universe changed.")
}

released_states <- read_hdf5_category(input, "cell_type")[cell_index]
if (anyNA(released_states) || any(released_states != current$state)) {
  stop("Released cell states differ from the current cell manifest.")
}
metadata <- data.frame(
  state = released_states,
  donor_id = as.character(h5read(input, "/obs/DonorID")[cell_index]),
  site = read_hdf5_category(input, "Site")[cell_index],
  batch = read_hdf5_category(input, "batch")[cell_index],
  sample = read_hdf5_category(input, "Samplename")[cell_index],
  row.names = current$cell_id,
  check.names = FALSE
)

message("Loading the complete current erythroid GEX count matrix from H5AD...")
h5_counts <- H5SparseMatrix(input, "/layers/counts")
counts <- as(h5_counts[gex_index, cell_index], "dgCMatrix")
rownames(counts) <- gex_names
colnames(counts) <- current$cell_id
if (!identical(dim(counts), c(13431L, 10952L)) || any(counts@x < 0)) {
  stop("Extracted erythroid GEX count matrix is invalid.")
}

message("Recomputing the erythroid-only GEX manifold...")
object <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata,
  min.cells = 0,
  min.features = 0,
  project = "GSE194122_erythroid"
)
object <- NormalizeData(object, normalization.method = "LogNormalize", verbose = FALSE)
object <- FindVariableFeatures(
  object,
  selection.method = "vst", nfeatures = 3000L, verbose = FALSE
)
object <- ScaleData(object, features = VariableFeatures(object), verbose = FALSE)
object <- RunPCA(
  object,
  features = VariableFeatures(object), npcs = 50L,
  seed.use = 20260822L, verbose = FALSE
)
object <- harmony::RunHarmony(
  object,
  group.by.vars = "batch",
  reduction.use = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  plot_convergence = FALSE,
  verbose = TRUE
)
object <- FindNeighbors(
  object,
  reduction = "harmony", dims = 1:30,
  k.param = 30L, verbose = FALSE
)
object <- RunUMAP(
  object,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 30L,
  min.dist = 0.3,
  metric = "cosine",
  seed.use = 20260822L,
  reduction.name = "erythroid_umap",
  reduction.key = "ERUMAP_",
  verbose = TRUE
)

embedding <- as.data.table(
  Embeddings(object, "erythroid_umap"),
  keep.rownames = "cell_id"
)
setnames(embedding, c("ERUMAP_1", "ERUMAP_2"), c("UMAP_1", "UMAP_2"))
embedding <- merge(
  embedding,
  as.data.table(metadata, keep.rownames = "cell_id"),
  by = "cell_id", all.x = TRUE, sort = FALSE
)
if (nrow(embedding) != 10952L || any(!is.finite(embedding$UMAP_1)) ||
  any(!is.finite(embedding$UMAP_2)) || anyNA(embedding$state)) {
  stop("The recomputed erythroid UMAP failed its completeness checks.")
}
setorder(embedding, cell_id)

fwrite(
  embedding,
  file.path(stage, "erythroid_gex_umap.tsv.gz"),
  sep = "\t", quote = FALSE
)
fwrite(
  embedding[, .(
    cells = .N,
    donors = uniqueN(donor_id),
    batches = uniqueN(batch),
    UMAP_1_median = median(UMAP_1),
    UMAP_2_median = median(UMAP_2)
  ), by = state],
  file.path(stage, "state_embedding_summary.tsv"),
  sep = "\t", quote = FALSE
)
fwrite(
  data.table(
    component = seq_along(Stdev(object, reduction = "pca")),
    standard_deviation = Stdev(object, reduction = "pca")
  ),
  file.path(stage, "pca_standard_deviation.tsv"),
  sep = "\t", quote = FALSE
)
if (!file.rename(stage, output)) stop("Unable to promote erythroid GEX UMAP outputs.")
cat("gse194122_erythroid_gex_umap\tPASS\t", output, "\n", sep = "")
