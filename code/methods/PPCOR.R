if (!requireNamespace("ppcor", quietly = TRUE)) {
  pak::pak("ppcor")
}

run_PPCOR <- function(
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
  results <- try(
    ppcor::pcor(
      x = t(matrix),
      method = "spearman"
    )
  )

  if (class(results) == "try-error") {
    return(NULL)
  }
  estimate <- as.matrix(results$estimate)
  edge_index <- which(row(estimate) != col(estimate), arr.ind = TRUE)
  network_table <- data.frame(
    regulator = gene_names[edge_index[, "row"]],
    target = gene_names[edge_index[, "col"]],
    weight = estimate[edge_index],
    stringsAsFactors = FALSE
  )
  network_table <- network_table[
    order(abs(network_table$weight), decreasing = TRUE), ,
    drop = FALSE
  ]
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
