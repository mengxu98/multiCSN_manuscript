if (!requireNamespace("SCORPION", quietly = TRUE)) {
  stop(
    "Package 'SCORPION' is required for the SCORPION benchmark. ",
    "Install it in the benchmark R environment before running this method."
  )
}

normalize_scorpion_matrix <- function(matrix) {
  matrix <- as.matrix(matrix)
  storage.mode(matrix) <- "numeric"

  if (is.null(rownames(matrix)) || any(!nzchar(rownames(matrix)))) {
    stop("SCORPION requires gene names as matrix rownames.")
  }
  if (is.null(colnames(matrix)) || any(!nzchar(colnames(matrix)))) {
    colnames(matrix) <- paste0("Cell", seq_len(ncol(matrix)))
  }

  matrix
}

# Method parameters are passed as function arguments.
scorpion_npc <- function(matrix, requested_npc = NULL) {
  if (is.null(requested_npc)) {
    requested_npc <- 25L
  }
  if (!is.finite(requested_npc) || requested_npc <= 0) {
    requested_npc <- 25L
  }

  max_npc <- max(1L, min(nrow(matrix), ncol(matrix)) - 1L)
  min(as.integer(requested_npc), max_npc)
}

scorpion_edge_table <- function(reg_net) {
  if (is.null(dim(reg_net))) {
    stop("SCORPION returned a regulation network without matrix dimensions.")
  }

  reg_net <- as.matrix(reg_net)
  regulator_names <- rownames(reg_net)
  target_names <- colnames(reg_net)
  if (is.null(regulator_names)) {
    regulator_names <- paste0("regulator_", seq_len(nrow(reg_net)))
  }
  if (is.null(target_names)) {
    target_names <- paste0("target_", seq_len(ncol(reg_net)))
  }

  index <- which(!is.na(reg_net), arr.ind = TRUE)
  network_table <- data.frame(
    regulator = regulator_names[index[, "row"]],
    target = target_names[index[, "col"]],
    weight = as.numeric(reg_net[index]),
    stringsAsFactors = FALSE
  )
  network_table <- network_table[network_table$regulator != network_table$target, , drop = FALSE]
  network_table <- network_table[order(network_table$weight, decreasing = TRUE), , drop = FALSE]
  rownames(network_table) <- NULL
  network_table
}

pad_scorpion_matrix <- function(matrix, min_genes = 10L) {
  if (nrow(matrix) >= min_genes) {
    return(list(matrix = matrix, original_genes = rownames(matrix)))
  }

  original_genes <- rownames(matrix)
  pad_count <- min_genes - nrow(matrix)
  cell_signal <- colMeans(matrix, na.rm = TRUE)
  if (anyNA(cell_signal)) {
    cell_signal[is.na(cell_signal)] <- 0
  }
  scale_value <- stats::sd(as.numeric(matrix), na.rm = TRUE)
  if (!is.finite(scale_value) || scale_value <= 0) {
    scale_value <- 1
  }

  pad_matrix <- vapply(seq_len(pad_count), function(idx) {
    cell_signal + (idx * seq_along(cell_signal) * scale_value * 1e-9)
  }, numeric(ncol(matrix)))
  pad_matrix <- t(pad_matrix)
  rownames(pad_matrix) <- paste0(".SCORPION_PAD_", seq_len(pad_count))
  colnames(pad_matrix) <- colnames(matrix)

  list(
    matrix = rbind(matrix, pad_matrix),
    original_genes = original_genes
  )
}

run_SCORPION <- function(
  matrix,
  regulators = NULL,
  targets = NULL,
  filename = "",
  n_cores = NULL,
  n_pc = NULL,
  gamma_value = NULL,
  show_progress = FALSE
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

  matrix <- normalize_scorpion_matrix(matrix)
  padded <- pad_scorpion_matrix(matrix)
  matrix <- padded$matrix

  if (is.null(n_cores)) {
    n_cores <- 1L
  }
  if (!is.finite(n_cores) || n_cores <= 0) {
    n_cores <- 1L
  }

  if (is.null(gamma_value)) {
    gamma_value <- 10
  }
  if (!is.finite(gamma_value) || gamma_value <= 0) {
    gamma_value <- 10
  }

  n_pc <- scorpion_npc(matrix, n_pc)

  result <- try(
    SCORPION::scorpion(
      tfMotifs = NULL,
      gexMatrix = matrix,
      ppiNet = NULL,
      nCores = as.integer(n_cores),
      nPC = as.integer(n_pc),
      gammaValue = gamma_value,
      showProgress = isTRUE(show_progress)
    ),
    silent = TRUE
  )

  if (inherits(result, "try-error")) {
    return(NULL)
  }

  reg_net <- if (is.list(result) && !is.null(result$regNet)) {
    result$regNet
  } else {
    result
  }
  network_table <- scorpion_edge_table(reg_net)
  network_table <- network_table[
    network_table$regulator %in% padded$original_genes &
      network_table$target %in% padded$original_genes, ,
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

  network_table
}
