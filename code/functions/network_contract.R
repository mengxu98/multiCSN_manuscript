source("code/functions/network_io.R")

# Validate actual endpoint tables, not a directory or completion marker.
validate_network_result <- function(root, expected = NULL) {
  endpoints <- c(`TF-gene` = "tf_gene.tsv", `TF-region` = "tf_region.tsv",
                 `region-gene` = "region_gene.tsv", `factorized-chain` = "triplets.tsv",
                 `mediated-TF-gene` = "mediated_tf_gene.tsv")
  required <- c(unname(endpoints), "endpoint_summary.tsv", "tf_gene_accounting.tsv",
                "tf_region_accounting.tsv", "region_gene_accounting.tsv")
  paths <- file.path(root, required)
  if (any(!file.exists(paths)) || any(file.info(paths)$size == 0))
    stop("Incomplete network output: ", root)
  summary <- data.table::fread(file.path(root, "endpoint_summary.tsv"))
  if (!all(c("endpoint", "selected") %in% names(summary)) ||
      anyDuplicated(summary$endpoint) || !setequal(summary$endpoint, names(endpoints)))
    stop("Invalid endpoint summary: ", root)
  for (endpoint in names(endpoints)) {
    tab <- data.table::fread(file.path(root, endpoints[[endpoint]]))
    keys <- switch(endpoint, `TF-gene` = c("regulator", "target"),
      `TF-region` = c("regulator", "region"), `region-gene` = c("region", "target"),
      `factorized-chain` = c("regulator", "region", "target"),
      `mediated-TF-gene` = c("regulator", "target"))
    if (!all(keys %in% names(tab)) || anyNA(tab[, ..keys]) || anyDuplicated(tab[, ..keys]))
      stop("Invalid endpoint keys: ", endpoint, " in ", root)
    if (nrow(tab) != summary$selected[match(endpoint, summary$endpoint)])
      stop("Endpoint count mismatch: ", endpoint, " in ", root)
  }
  if (!is.null(expected)) {
    provenance_path <- file.path(root, "run_provenance.rds")
    if (!file.exists(provenance_path))
      stop("No verifiable run identity; archive and rerun rather than assume PASS: ", root)
    provenance <- readRDS(provenance_path)
    for (key in names(expected)) {
      value <- expected[[key]]
      if (!is.null(value) && key %in% c("input", "validation", "cell_manifest"))
        value <- normalizePath(value, mustWork = TRUE)
      if (!isTRUE(all.equal(provenance$settings[[key]], value, check.attributes = FALSE)))
        stop("Run identity mismatch: ", key)
    }
    files <- provenance$inputs$file
    if (!length(files) || any(!file.exists(files)) ||
        !identical(unname(sha256_files(files)), provenance$inputs$sha256))
      stop("Input/code/package identity changed: ", root)
    files <- file.path(root, provenance$outputs$file)
    if (!all(required %in% provenance$outputs$file) || any(!file.exists(files)) ||
        !identical(unname(sha256_files(files)), provenance$outputs$sha256))
      stop("Output content changed: ", root)
  }
  invisible(TRUE)
}
