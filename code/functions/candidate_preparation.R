read_h5ad_csr_transposed <- function(path, group, features, cells) {
  attributes <- rhdf5::h5readAttributes(path, group)
  if (!identical(as.character(attributes[["encoding-type"]]), "csr_matrix") ||
      !identical(as.integer(attributes$shape), c(length(cells), length(features)))) {
    stop("Expected the released cell-by-feature CSR matrix: ", group)
  }
  methods::new("dgCMatrix", Dim = as.integer(c(length(features), length(cells))),
    Dimnames = list(features, cells),
    x = as.numeric(rhdf5::h5read(path, paste0(group, "/data"))),
    i = as.integer(rhdf5::h5read(path, paste0(group, "/indices"))),
    p = as.integer(rhdf5::h5read(path, paste0(group, "/indptr"))))
}

read_released_h5ad_column <- function(path, axis, field) {
  values <- as.vector(rhdf5::h5read(path, paste0("/", axis, "/", field)))
  category <- paste0("/", axis, "/__categories/", field)
  handle <- rhdf5::H5Fopen(path)
  has_category <- tryCatch(rhdf5::H5Lexists(handle, category), finally = rhdf5::H5Fclose(handle))
  if (has_category) {
    levels <- as.character(rhdf5::h5read(path, category))
    if (anyNA(values) || any(values < -1L | values >= length(levels))) stop("Invalid H5AD category: ", field)
    result <- rep(NA_character_, length(values))
    present <- values >= 0L
    result[present] <- levels[values[present] + 1L]
    return(result)
  }
  values
}

human_hg38_annotation <- function() {
  annotation <- Signac::GetGRangesFromEnsDb(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
  mapping <- stats::setNames(c(paste0("chr", 1:22), "chrX", "chrY", "chrM"), c(as.character(1:22), "X", "Y", "MT"))
  mapping <- mapping[names(mapping) %in% GenomeInfoDb::seqlevels(annotation)]
  annotation <- GenomeInfoDb::renameSeqlevels(annotation, mapping)
  GenomeInfoDb::genome(annotation) <- "hg38"
  annotation
}

create_full_chromatin_assay <- function(counts, annotation) {
  ranges <- multiCSN::parse_peak_ranges(rownames(counts))
  GenomeInfoDb::genome(ranges) <- "hg38"
  Signac::CreateChromatinAssay(counts = counts, ranges = ranges, annotation = annotation,
    genome = GenomeInfoDb::seqinfo(ranges), min.cells = 0, min.features = 0)
}

read_paired_sparse_h5 <- function(path, group) {
  methods::new("dgCMatrix", Dim = as.integer(rhdf5::h5read(path, paste0(group, "/shape"))),
    Dimnames = list(as.character(rhdf5::h5read(path, paste0(group, "/features"))),
                    as.character(rhdf5::h5read(path, paste0(group, "/barcodes")))),
    x = as.numeric(rhdf5::h5read(path, paste0(group, "/data"))),
    i = as.integer(rhdf5::h5read(path, paste0(group, "/indices"))),
    p = as.integer(rhdf5::h5read(path, paste0(group, "/indptr"))))
}

read_gse162170_release <- function(input, annotation = human_hg38_annotation()) {
  force(annotation)
  python <- Sys.which("python3")
  if (!nzchar(python)) stop("Python with pandas, scipy and h5py is required.")
  sparse_path <- tempfile("gse162170-paired-", fileext = ".h5")
  on.exit(unlink(sparse_path), add = TRUE)
  script <- normalizePath("code/functions/paired_counts.py", mustWork = TRUE)
  expression <- paste0("import runpy; m=runpy.run_path(", encodeString(script, quote = '"'),
    "); m['prepare_paired_counts'](", encodeString(normalizePath(input, mustWork = TRUE), quote = '"'),
    ", ", encodeString(sparse_path, quote = '"'), ")")
  status <- system2(python, c("-B", "-c", shQuote(expression)))
  if (status != 0L) stop("Unable to read the released paired count matrices.")
  build_gse162170_object(sparse_path, input, annotation)
}

build_gse162170_object <- function(sparse_path, input, annotation) {
  rna <- read_paired_sparse_h5(sparse_path, "rna")
  symbols <- as.character(AnnotationDbi::mapIds(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
    keys = rownames(rna), keytype = "GENEID", column = "GENENAME", multiVals = "first"))
  missing <- is.na(symbols) | !nzchar(symbols)
  symbols[missing] <- rownames(rna)[missing]
  symbols <- gsub("_", "-", symbols, fixed = TRUE)
  unique_symbols <- unique(symbols)
  mapping <- Matrix::sparseMatrix(i = match(symbols, unique_symbols), j = seq_along(symbols),
    x = 1, dims = c(length(unique_symbols), length(symbols)), dimnames = list(unique_symbols, rownames(rna)))
  rna <- mapping %*% rna
  metadata <- as.data.frame(data.table::fread(file.path(input, "GSE162170_multiome_cell_metadata.txt.gz")))
  rownames(metadata) <- metadata$Cell.ID
  if (!identical(colnames(rna), rownames(metadata))) stop("RNA cells do not match the paired metadata.")
  clusters <- data.table::fread(file.path(input, "GSE162170_multiome_cluster_names.txt.gz"))
  clusters <- clusters[Assay == "Multiome RNA"]
  metadata$author_cluster_name <- clusters$Cluster.Name[match(metadata$seurat_clusters, clusters$Cluster.ID)]
  if (anyNA(metadata$author_cluster_name)) stop("Paired author cluster names do not map exactly.")
  metadata$same_cell_measurement <- TRUE
  metadata$alignment_unit <- metadata$Cell.ID
  object <- SeuratObject::CreateSeuratObject(rna, project = "GSE162170_full", meta.data = metadata, min.cells = 0, min.features = 0)
  atac <- read_paired_sparse_h5(sparse_path, "atac")
  if (!identical(colnames(atac), colnames(object))) stop("ATAC cells do not match the paired RNA cells.")
  object[["ATAC"]] <- create_full_chromatin_assay(atac, annotation)
  object <- Seurat::NormalizeData(object, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  object <- SeuratObject::SetAssayData(object, assay = "ATAC", layer = "data",
    new.data = Signac::RunTFIDF(SeuratObject::LayerData(object, assay = "ATAC", layer = "counts"), method = 1, scale.factor = 10000, verbose = FALSE))
  object$nCount_peaks <- object$nFeature_ATAC
  SeuratObject::DefaultAssay(object) <- "ATAC"
  object
}

read_gse194122_release <- function(input, annotation = human_hg38_annotation()) {
  cells <- as.character(rhdf5::h5read(input, "/obs/_index"))
  features <- as.character(rhdf5::h5read(input, "/var/_index"))
  feature_types <- read_released_h5ad_column(input, "var", "feature_types")
  if (!length(cells) || anyDuplicated(cells) || anyDuplicated(features) ||
      length(feature_types) != length(features) || anyNA(feature_types) ||
      !all(c("GEX", "ATAC") %in% feature_types)) {
    stop("Invalid released cell or feature identifiers.")
  }
  obs <- rhdf5::h5ls(input, recursive = TRUE)
  fields <- obs$name[obs$group == "/obs" & obs$otype == "H5I_DATASET" & obs$name != "_index"]
  metadata <- as.data.frame(stats::setNames(lapply(fields, function(field) read_released_h5ad_column(input, "obs", field)), fields), stringsAsFactors = FALSE)
  rownames(metadata) <- cells
  metadata$cell_id <- cells
  counts <- read_h5ad_csr_transposed(input, "/layers/counts", features, cells)
  object <- SeuratObject::CreateSeuratObject(counts = counts[feature_types == "GEX", , drop = FALSE],
    project = "GSE194122_full", meta.data = metadata, min.cells = 0, min.features = 0)
  object[["ATAC"]] <- create_full_chromatin_assay(counts[feature_types == "ATAC", , drop = FALSE], annotation)
  rm(counts)
  object <- SeuratObject::SetAssayData(object, assay = "ATAC", layer = "data",
    new.data = Signac::RunTFIDF(SeuratObject::LayerData(object, assay = "ATAC", layer = "counts"),
                              method = 1, scale.factor = 10000, verbose = FALSE))
  object <- Seurat::NormalizeData(object, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  object$nCount_peaks <- object$nFeature_ATAC
  SeuratObject::DefaultAssay(object) <- "ATAC"
  object
}

align_candidate_annotation <- function(object) {
  annotation <- Signac::Annotation(object[["ATAC"]])
  chromosomes <- GenomeInfoDb::seqlevels(multiCSN::parse_peak_ranges(rownames(object[["ATAC"]])))
  GenomeInfoDb::seqlevels(annotation) <- union(GenomeInfoDb::seqlevels(annotation), chromosomes)
  GenomeInfoDb::genome(annotation) <- "hg38"
  Signac::Annotation(object[["ATAC"]]) <- annotation
  object
}

prepare_dataset_candidates <- function(dataset, input, destination,
    motif_dir = "data/resources/motifs") {
  suppressPackageStartupMessages({library(Matrix); library(SeuratObject); library(Signac); library(multiCSN); library(motifmatchr)})
  if (!dataset %in% c("GSE162170", "GSE194122")) stop("Unknown candidate dataset.")
  if (file.exists(destination)) stop("Candidate output already exists: ", destination)
  object <- if (dataset == "GSE162170") read_gse162170_release(input) else read_gse194122_release(input)
  rna <- rownames(object[["RNA"]])
  object <- align_candidate_annotation(object)
  object <- multiCSN::initiate_object(object, group.by = "aggregate", targets = rna,
    filter_mode = "unfiltered", filter_by = "aggregate", rna_assay = "RNA", peak_assay = "ATAC",
    regions = NULL, regions_extend = 0, exclude_exons = TRUE, verbose = TRUE)
  load_one <- function(path) {
    environment <- new.env(parent = emptyenv()); names <- load(path, envir = environment)
    if (length(names) != 1L) stop("Expected one motif resource object.")
    environment[[names[[1L]]]]
  }
  object <- multiCSN::find_motifs(object, pfm = load_one(file.path(motif_dir, "motifs.rda")),
    motif_tfs = load_one(file.path(motif_dir, "motif2tf.rda")),
    genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, backend = "motifmatchr", verbose = TRUE)
  selected <- multiCSN::get_attribute(object, celltypes = "aggregate", attribute = "peaks")
  genes <- multiCSN::get_attribute(object, celltypes = "aggregate", attribute = "genes")
  motifs <- multiCSN::NetworkRegions(object)@motifs@data
  if (!identical(sort(genes), sort(rna)) || !setequal(rownames(motifs), selected))
    stop("Candidate gene or peak names do not match the input.")
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(".candidates-", tmpdir = dirname(destination))
  on.exit(unlink(stage), add = TRUE)
  saveRDS(object, stage, compress = FALSE)
  if (file.exists(destination)) stop("Candidate output appeared during preparation and will not be overwritten.")
  if (!file.rename(stage, destination)) stop("Unable to publish the candidate object.")
  invisible(destination)
}
