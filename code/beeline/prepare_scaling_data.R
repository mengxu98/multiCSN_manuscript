build_scaling_case <- function(total_genes, tf_count, cell_count, seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  pseudotime <- seq(0, 1, length.out = cell_count)
  basis <- cbind(pseudotime, 1 - pseudotime,
    exp(-((pseudotime - 0.5) / 0.18)^2),
    1 / (1 + exp(-10 * (pseudotime - 0.55))),
    0.5 + 0.5 * sin(2 * pi * pseudotime))
  regulators <- matrix(0, tf_count, cell_count)
  for (i in seq_len(tf_count)) {
    regulators[i, ] <- pmax(0.2 + (i %% 11) * 0.03 +
      (0.8 + (i %% 7) * 0.12) * basis[, ((i - 1L) %% ncol(basis)) + 1L] +
      rnorm(cell_count, 0, 0.08 + (i %% 5) * 0.01), 0)
  }
  target_count <- total_genes - tf_count
  targets <- matrix(0, target_count, cell_count)
  edges <- vector("list", target_count)
  for (i in seq_len(target_count)) {
    size <- sample(c(1L, 2L), 1L, prob = c(0.7, 0.3))
    indices <- sort(sample.int(tf_count, size))
    weights <- runif(size, 0.45, 1.1)
    directions <- sample(c(-1, 1), size, replace = TRUE, prob = c(0.25, 0.75))
    signal <- colSums(regulators[indices, , drop = FALSE] * (weights * directions))
    trend <- (0.08 + (i %% 9) * 0.015) * pseudotime^((i %% 4) + 1L)
    targets[i, ] <- pmax(0.15 + trend + signal + rnorm(cell_count, 0, 0.15 + (i %% 6) * 0.01), 0)
    edges[[i]] <- data.frame(Gene1 = sprintf("TF%04d", indices), Gene2 = sprintf("G%04d", i),
                            Type = ifelse(directions > 0, "+", "-"))
  }
  expression <- rbind(regulators, targets)
  rownames(expression) <- c(sprintf("TF%04d", seq_len(tf_count)), sprintf("G%04d", seq_len(target_count)))
  colnames(expression) <- sprintf("Cell%05d", seq_len(cell_count))
  list(expression = expression,
       pseudotime = data.frame(PseudoTime = pseudotime, row.names = colnames(expression)),
       reference = do.call(rbind, edges))
}

scaling_design <- function() {
  cells <- c(100L, 500L, 1000L, 2000L, 5000L)
  targets <- c(50L, 100L, 200L, 500L, 1000L)
  rbind(data.frame(Case = sprintf("scaling_cells_%d_genes_1000_tfs_500", cells),
    Scenario = "cell_scaling", RequestedCells = cells, RequestedGenes = 1000L,
    RequestedTFs = 500L, TargetGeneCount = 500L, TotalGeneCount = 1000L),
    data.frame(Case = sprintf("scaling_genes_%d_cells_5000_tfs_500", targets),
    Scenario = "gene_scaling", RequestedCells = 5000L, RequestedGenes = targets,
    RequestedTFs = 500L, TargetGeneCount = targets, TotalGeneCount = targets + 500L))
}

prepare_scaling_data <- function(output = "data/beeline/scaling") {
  design <- scaling_design()
  if (dir.exists(output)) {
    manifest <- file.path(output, "scaling_manifest.csv")
    if (!identical(read.csv(manifest), design)) stop("Scaling design does not match.")
    for (i in seq_len(nrow(design))) {
      required <- file.path(output, design$Scenario[i], design$Case[i],
                            c("ExpressionData.csv", "PseudoTime.csv", "refNetwork.csv", "case_metadata.csv"))
      if (any(!file.exists(required))) stop("Existing scaling input is incomplete: ", design$Case[i])
    }
    return(invisible(output))
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(".scaling-inputs-", tmpdir = dirname(output))
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  for (i in seq_len(nrow(design))) {
    row <- design[i, , drop = FALSE]
    seed <- if (row$Scenario == "cell_scaling") 1000L + i else 2000L + i - 5L
    case <- build_scaling_case(row$TotalGeneCount, row$RequestedTFs, row$RequestedCells, seed)
    path <- file.path(stage, row$Scenario, row$Case)
    dir.create(path, recursive = TRUE)
    write.csv(case$expression, file.path(path, "ExpressionData.csv"))
    write.csv(case$pseudotime, file.path(path, "PseudoTime.csv"))
    write.csv(case$reference, file.path(path, "refNetwork.csv"), row.names = FALSE)
    write.csv(row, file.path(path, "case_metadata.csv"), row.names = FALSE)
  }
  write.csv(design, file.path(stage, "scaling_manifest.csv"), row.names = FALSE)
  if (!file.rename(stage, output)) stop("Unable to publish the scaling inputs.")
  invisible(output)
}

if (sys.nframe() == 0L) prepare_scaling_data()
