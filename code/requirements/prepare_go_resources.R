prepare_go_resources <- function(resource_root = "data/resources/fig5_go_3.22.0") {
  for (package in c("GO.db", "org.Hs.eg.db")) {
    if (!requireNamespace(package, quietly = TRUE) || as.character(utils::packageVersion(package)) != "3.22.0") {
      stop("Install ", package, " 3.22.0 from Bioconductor 3.22 before preparing the fixed GO resource.")
    }
  }
  cache <- file.path(resource_root, "cache")
  path <- file.path(cache, "af3b10e4c2379820a80657c5a6685dc4.Rcache")
  if (file.exists(path)) {
    resource <- R.cache::loadCache(pathname = path)
    if (as.character(resource$version) != "3.22.0" ||
        !all(c("Term", "symbol") %in% names(resource$TERM2GENE)) ||
        !all(c("Term", "Name", "ONTOLOGY") %in% names(resource$TERM2NAME)) ||
        !nrow(resource$TERM2GENE) || !nrow(resource$TERM2NAME) ||
        !identical(unique(as.character(resource$TERM2NAME$ONTOLOGY)), "BP")) {
      stop("The fixed GO cache failed validation.")
    }
    return(invisible(path))
  }
  dir.create(resource_root, recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(".go-resource-", tmpdir = resource_root)
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  previous <- R.cache::getCacheRootPath()
  on.exit(R.cache::setCacheRootPath(previous), add = TRUE)
  R.cache::setCacheRootPath(stage)
  resource <- scop::PrepareDB(species = "Homo_sapiens", db = "GO_BP", db_update = TRUE,
                              db_IDtypes = "symbol", verbose = FALSE)[["Homo_sapiens"]][["GO_BP"]]
  R.cache::saveCache(resource, key = list(resource$version, "Homo_sapiens", "GO_BP"))
  built <- file.path(stage, basename(path))
  if (!file.exists(built)) stop("GO cache key differs from the declared resource.")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(built, path)) stop("Unable to publish the GO annotation resource.")
  invisible(path)
}

prepare_display_go_similarity <- function(resource_root = "data/resources/fig5_go_3.22.0",
    enrichment_path = "results/gse194122/donor_stability/stable_rewiring_GO_BP.tsv.gz") {
  path <- prepare_go_resources(resource_root)
  resource <- R.cache::loadCache(pathname = path)
  x <- data.table::fread(enrichment_path)
  a <- x[regulator == "HMGA1" & transition_order == 4 & topology_class == "lost" & p_adjust < 0.05, go_id]
  b <- intersect(x[regulator == "BNC2" & transition_order == 3 & topology_class == "gained" & p_adjust < 0.05, go_id],
                 x[regulator == "BNC2" & transition_order == 4 & topology_class == "lost" & p_adjust < 0.05, go_id])
  terms <- sort(union(a, b))
  if (!length(terms) || anyNA(terms) || anyDuplicated(terms)) stop("Invalid display GO universe.")
  destination <- file.path(resource_root, "display_wang_similarity.rds")
  if (file.exists(destination)) {
    saved <- readRDS(destination)
    if (!setequal(rownames(saved), terms) || !setequal(colnames(saved), terms)) {
      stop("Saved GO similarity terms do not match the enrichment input.")
    }
    if (!is.numeric(saved) || any(!is.finite(saved)) || any(saved < 0 | saved > 1)) {
      stop("Invalid GO similarities.")
    }
    return(invisible(destination))
  }
  similarity <- GOSemSim::mgoSim(terms, terms, semData = resource$semData, measure = "Wang", combine = NULL)
  stage <- tempfile(".wang-", tmpdir = resource_root)
  on.exit(unlink(stage), add = TRUE)
  saveRDS(similarity, stage)
  if (!file.rename(stage, destination)) stop("Unable to publish the Wang similarity resource.")
  invisible(destination)
}

if (sys.nframe() == 0L) prepare_go_resources()
