source("code/functions/network_io.R")

# Validate actual endpoint tables, not a directory or completion marker.
validate_network_result <- function(root) {
  spec <- multiCSN::network_endpoints()
  endpoints <- stats::setNames(spec$file, spec$endpoint)
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
    tab <- multiCSN::read_network_endpoints(root, endpoints = endpoint)[[endpoint]]
    if (nrow(tab) != summary$selected[match(endpoint, summary$endpoint)])
      stop("Endpoint count mismatch: ", endpoint, " in ", root)
  }
  invisible(TRUE)
}
