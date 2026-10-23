source("code/functions/network_io.R")

# Validate actual endpoint tables, not a directory or completion marker.
validate_network_result <- function(root) {
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
  invisible(TRUE)
}
