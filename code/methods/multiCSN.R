csn_export <- function(name) {
  if (!requireNamespace("inferCSN", quietly = TRUE)) {
    stop("Install inferCSN before running the benchmark.", call. = FALSE)
  }
  getExportedValue("inferCSN", name)
}

run_multiCSN <- function(
  matrix,
  pseudotime,
  regulators = NULL,
  targets = NULL,
  max_support_size = NULL,
  lag_fraction = 0.05,
  lag_steps = NULL,
  cores = 1,
  filename = "",
  # Default FALSE: verbose logging in inferCSN incurs overhead (~20 ms per call)
  # which skews benchmarking on small test cases. Pass verbose = TRUE if progress is desired.
  verbose = FALSE
) {
  csn_fn <- csn_export("inferCSN")

  network_table <- try(
    csn_fn(
      matrix,
      pseudotime = pseudotime,
      regulators = regulators,
      targets = targets,
      max_support_size = max_support_size,
      lag_fraction = lag_fraction,
      lag_steps = lag_steps,
      cores = cores,
      verbose = verbose
    )
  )
  if (inherits(network_table, "try-error")) {
    return(NULL)
  }

  if (nzchar(filename)) {
    write.csv(network_table, file = filename, row.names = FALSE)
  }
  network_table
}
