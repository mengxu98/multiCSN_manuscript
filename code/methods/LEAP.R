if (!requireNamespace("LEAP", quietly = TRUE)) {
  stop(
    "Package 'LEAP' is required for the LEAP benchmark. ",
    "Install it in the benchmark R environment before running this method."
  )
}

run_LEAP <- function(
  matrix,
  pseudotime,
  regulators = NULL,
  targets = NULL,
  filename = ""
) {
  matrix <- as.matrix(matrix)
  if (is.null(rownames(matrix)) || is.null(colnames(matrix)) ||
      anyDuplicated(rownames(matrix)) || anyDuplicated(colnames(matrix))) {
    stop("LEAP requires unique gene and cell identifiers.")
  }
  if (is.null(dim(pseudotime))) {
    pseudotime <- matrix(pseudotime, ncol = 1L,
                         dimnames = list(names(pseudotime), "trajectory"))
  } else {
    pseudotime <- as.matrix(pseudotime)
  }
  if (!is.numeric(pseudotime) || nrow(pseudotime) != ncol(matrix)) {
    stop("LEAP pseudotime must be numeric with one row per expression cell.")
  }
  if (!is.null(rownames(pseudotime))) {
    if (anyDuplicated(rownames(pseudotime)) ||
        !setequal(rownames(pseudotime), colnames(matrix))) {
      stop("LEAP pseudotime and expression cell identifiers do not match.")
    }
    pseudotime <- pseudotime[colnames(matrix), , drop = FALSE]
  }
  gene_keep <- rownames(matrix)
  if (!is.null(regulators) && length(regulators) > 0) {
    gene_keep <- intersect(gene_keep, regulators)
  }
  if (!is.null(targets) && length(targets) > 0) {
    gene_keep <- union(gene_keep, intersect(rownames(matrix), targets))
  }
  matrix <- matrix[unique(gene_keep), , drop = FALSE]
  if (nrow(matrix) < 2L) stop("LEAP requires at least two candidate genes.")
  branches <- lapply(seq_len(ncol(pseudotime)), function(branch) {
    cells <- which(is.finite(pseudotime[, branch]))
    cells[order(pseudotime[cells, branch], colnames(matrix)[cells], method = "radix")]
  })
  branches <- unique(branches[lengths(branches) >= 3L])
  if (!length(branches)) stop("LEAP requires a trajectory with at least three finite-time cells.")
  networks <- lapply(branches, function(cells) {
    results <- LEAP::MAC_counter(
      data = matrix[, cells, drop = FALSE], MAC_cutoff = 0,
      lag_matrix = FALSE, symmetric = FALSE
    )
    data.frame(regulator = rownames(matrix)[results[, "Row gene index"]],
               target = rownames(matrix)[results[, "Column gene index"]],
               weight = as.numeric(results[, "Correlation"]))
  })
  network_table <- do.call(rbind, networks)
  network_table <- network_table[is.finite(network_table$weight), , drop = FALSE]
  network_table <- network_table[order(-abs(network_table$weight), method = "radix"), , drop = FALSE]
  network_table <- network_table[!duplicated(network_table[c("regulator", "target")]), , drop = FALSE]
  if (!is.null(regulators) && length(regulators) > 0) {
    network_table <- network_table[network_table$regulator %in% regulators, , drop = FALSE]
  }
  if (!is.null(targets) && length(targets) > 0) {
    network_table <- network_table[network_table$target %in% targets, , drop = FALSE]
  }

  if (nzchar(filename)) write.csv(network_table, file = filename, row.names = FALSE)

  return(network_table)
}
