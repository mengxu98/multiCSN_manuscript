if (!requireNamespace("GENIE3", quietly = TRUE)) {
  pak::pak("GENIE3")
}

run_GENIE3 <- function(
  matrix,
  cores = 1,
  regulators = NULL,
  targets = NULL,
  filename = ""
) {
  network_table <- try(
    GENIE3::getLinkList(
      GENIE3::GENIE3(
        exprMatrix = matrix,
        regulators = regulators,
        targets = targets,
        nCores = cores
      )
    )
  )

  if (class(network_table) == "try-error") {
    return(NULL)
  }

  if (!is.null(regulators) && length(regulators) > 0) {
    network_table <- network_table[network_table[[1]] %in% regulators, , drop = FALSE]
  }
  if (!is.null(targets) && length(targets) > 0) {
    network_table <- network_table[network_table[[2]] %in% targets, , drop = FALSE]
  }

  write.csv(
    network_table,
    file = filename,
    row.names = FALSE
  )

  return(network_table)
}
