if (!requireNamespace("LEAP", quietly = TRUE)) {
  stop(
    "Package 'LEAP' is required for the LEAP benchmark. ",
    "Install it in the benchmark R environment before running this method."
  )
}

run_LEAP <- function(
  matrix,
  regulators = NULL,
  targets = NULL,
  filename = ""
) {
  gene_keep <- rownames(matrix)
  if (!is.null(regulators) && length(regulators) > 0) {
    gene_keep <- intersect(gene_keep, regulators)
  }
  if (!is.null(targets) && length(targets) > 0) {
    gene_keep <- union(gene_keep, intersect(rownames(matrix), targets))
  }
  matrix <- matrix[unique(gene_keep), , drop = FALSE]
  if (nrow(matrix) < 2) {
    return(NULL)
  }

  gene_names <- rownames(matrix)
  rownames(matrix) <- c()
  results <- try(
    LEAP::MAC_counter(
      data = matrix,
      MAC_cutoff = 0,
      lag_matrix = FALSE,
      symmetric = FALSE
    )
  )
  if (any(class(results) == "try-error")) {
    return(NULL)
  }
  gene1 <- gene_names[results[, "Row gene index"]]
  gene2 <- gene_names[results[, "Column gene index"]]
  score <- results[, "Correlation"]
  network_table <- data.frame(regulator = gene1, target = gene2, weight = score)
  if (!is.null(regulators) && length(regulators) > 0) {
    network_table <- network_table[network_table$regulator %in% regulators, , drop = FALSE]
  }
  if (!is.null(targets) && length(targets) > 0) {
    network_table <- network_table[network_table$target %in% targets, , drop = FALSE]
  }

  write.csv(
    network_table,
    file = filename,
    row.names = FALSE
  )

  return(network_table)
}
