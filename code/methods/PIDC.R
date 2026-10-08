ensure_pidc_available <- function() {
  if (requireNamespace("PIDC", quietly = TRUE)) {
    return("package")
  }

  vendor_candidates <- unique(c(
    file.path("vendor", "PIDC", "PIDC.R"),
    file.path("code", "vendor", "PIDC", "PIDC.R"),
    if (exists("benchmark_root_dir", mode = "function")) {
      file.path(benchmark_root_dir(), "vendor", "PIDC", "PIDC.R")
    } else {
      character(0)
    }
  ))
  for (vendor_path in vendor_candidates) {
    if (!file.exists(vendor_path)) {
      next
    }
    source(vendor_path, local = parent.frame())
    if (exists("PIDC", mode = "function", envir = parent.frame(), inherits = FALSE)) {
      return("vendor")
    }
  }

  stop("PIDC implementation not available via installed package or local vendor source.")

  invisible(NULL)
}

run_PIDC <- function(
  matrix,
  regulators = NULL,
  targets = NULL,
  cores = 1,
  filename = "",
  log_scale = FALSE,
  verbose = FALSE
) {
  pidc_backend <- ensure_pidc_available()
  matrix <- as.matrix(matrix)
  nonconstant_mask <- apply(matrix, 1, function(values) {
    values <- suppressWarnings(as.numeric(values))
    keep <- is.finite(values)
    if (sum(keep) < 2) {
      return(FALSE)
    }
    stats::sd(values[keep]) > 0
  })
  matrix <- matrix[nonconstant_mask, , drop = FALSE]
  if (nrow(matrix) < 2) {
    return(NULL)
  }

  weight_matrix <- try(
    if (identical(pidc_backend, "package")) {
      PIDC::PIDC(
        expMat = matrix,
        regulators = regulators,
        targets = targets,
        logScale = log_scale,
        ncores = cores,
        diag = "zero",
        verbose = verbose
      )
    } else {
      PIDC(
        expMat = matrix,
        regulators = regulators,
        targets = targets,
        logScale = log_scale,
        ncores = cores,
        diag = "zero",
        verbose = verbose
      )
    }
  )

  if (inherits(weight_matrix, "try-error")) {
    return(NULL)
  }

  network_table <- square_matrix_to_network(weight_matrix)
  if (nrow(network_table) == 0) {
    return(NULL)
  }
  if (!is.null(regulators) && length(regulators) > 0) {
    network_table <- network_table[
      network_table$regulator %in% regulators, ,
      drop = FALSE
    ]
  }
  if (!is.null(targets) && length(targets) > 0) {
    network_table <- network_table[
      network_table$target %in% targets, ,
      drop = FALSE
    ]
  }
  if (nrow(network_table) == 0) {
    return(NULL)
  }

  write.csv(
    network_table,
    file = filename,
    row.names = FALSE
  )

  network_table
}
