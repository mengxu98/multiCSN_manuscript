source("code/functions/runtime.R")

normalize_pseudotime <- function(pseudotime) {
  if (is.data.frame(pseudotime) || is.matrix(pseudotime)) {
    pseudotime <- as.data.frame(pseudotime, check.names = FALSE)
    pseudotime <- apply(pseudotime, 1, function(row_values) {
      numeric_values <- suppressWarnings(as.numeric(row_values))
      numeric_values <- numeric_values[!is.na(numeric_values)]
      if (length(numeric_values)) {
        numeric_values[[1]]
      } else {
        NA_real_
      }
    })
  }
  pseudotime <- as.numeric(pseudotime)
  if (anyNA(pseudotime)) {
    fallback <- seq(0, 1, length.out = length(pseudotime))
    pseudotime[is.na(pseudotime)] <- fallback[is.na(pseudotime)]
  }
  pseudotime <- pseudotime - min(pseudotime, na.rm = TRUE)
  max_time <- max(pseudotime, na.rm = TRUE)

  if (!is.finite(max_time) || max_time <= 0) {
    return(rep(0, length(pseudotime)))
  }

  pseudotime / max_time
}

find_scode_script <- function() {
  candidates <- c(
    file.path(benchmark_root_dir(), "vendor", "SCODE", "SCODE.R"),
    file.path(benchmark_root_dir(), "vendor", "SCODE", "scode.R")
  )
  candidates <- candidates[file.exists(candidates)]

  if (length(candidates) == 0) {
    stop(
      "SCODE script not found. Please place the official 'SCODE.R' script ",
      "under code/vendor/SCODE/."
    )
  }

  candidates[[1]]
}

run_SCODE <- function(
  matrix,
  pseudotime,
  regulators = NULL,
  targets = NULL,
  filename = "",
  latent_dim = NULL,
  n_iter = 100,
  n_runs = 5,
  rscript_bin = Sys.which("Rscript")
) {
  if (!nzchar(rscript_bin)) {
    stop("Rscript not found for running SCODE.")
  }

  script_path <- find_scode_script()

  gene_names <- rownames(matrix)
  matrix_input <- as.matrix(matrix)

  if (!is.null(regulators) && length(regulators) > 0) {
    keep_genes <- intersect(gene_names, regulators)
  } else {
    keep_genes <- gene_names
  }
  if (!is.null(targets) && length(targets) > 0) {
    keep_genes <- union(keep_genes, intersect(gene_names, targets))
  }
  keep_genes <- unique(keep_genes)

  matrix_input <- matrix_input[keep_genes, , drop = FALSE]
  if (nrow(matrix_input) < 2 || ncol(matrix_input) < 3) {
    return(NULL)
  }

  if (is.null(latent_dim)) {
    latent_dim <- min(4, nrow(matrix_input))
  }
  latent_dim <- max(1, min(latent_dim, nrow(matrix_input)))

  temp_dir <- tempfile("scode_")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  expression_file <- file.path(temp_dir, "expression.tsv")
  pseudotime_file <- file.path(temp_dir, "pseudotime.tsv")

  write.table(
    matrix_input,
    file = expression_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  pseudotime_df <- data.frame(
    cell = seq_len(ncol(matrix_input)) - 1,
    pseudotime = normalize_pseudotime(pseudotime)
  )
  write.table(
    pseudotime_df,
    file = pseudotime_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  network_mats <- vector("list", n_runs)

  for (run_id in seq_len(n_runs)) {
    output_dir <- file.path(temp_dir, paste0("run_", run_id))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    settings <- list(fdata = expression_file, ftime = pseudotime_file, dir = output_dir,
                     tfnum = nrow(matrix_input), pnum = latent_dim,
                     cnum = ncol(matrix_input), maxite = n_iter)
    expression <- r_function_expression(script_path, "run_scode", settings)

    run_log <- try(
      system2(
        command = rscript_bin,
        args = c("-e", shQuote(expression)),
        stdout = TRUE,
        stderr = TRUE
      ),
      silent = TRUE
    )

    if (inherits(run_log, "try-error")) {
      return(NULL)
    }
    status <- attr(run_log, "status")
    if (!is.null(status) && status != 0) {
      warning(
        "SCODE subprocess failed with status ", status, ".\n",
        paste(run_log, collapse = "\n")
      )
      return(NULL)
    }

    a_file <- file.path(output_dir, "A.txt")
    if (!file.exists(a_file)) {
      warning("SCODE did not produce expected file: ", a_file)
      return(NULL)
    }

    network_mat <- as.matrix(
      read.table(
        a_file,
        sep = "\t",
        header = FALSE,
        stringsAsFactors = FALSE
      )
    )
    if (!all(dim(network_mat) == c(nrow(matrix_input), nrow(matrix_input)))) {
      return(NULL)
    }

    rownames(network_mat) <- rownames(matrix_input)
    colnames(network_mat) <- rownames(matrix_input)
    network_mats[[run_id]] <- network_mat
  }

  averaged_mat <- Reduce(`+`, network_mats) / length(network_mats)
  network_table <- square_matrix_to_network(averaged_mat)

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
