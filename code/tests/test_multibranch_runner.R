source("code/methods/multiCSN.R")

captured_pseudotime <- NULL
csn_export <- function(name) {
  stopifnot(identical(name, "inferCSN"))
  function(matrix, pseudotime, ...) {
    captured_pseudotime <<- pseudotime
    data.frame(
      regulator = "tf",
      target = "target",
      weight = 1
    )
  }
}

pseudotime <- cbind(
  branch_a = c(1, 2, NA, NA),
  branch_b = c(NA, NA, 1, 2)
)
output_file <- tempfile(fileext = ".csv")
on.exit(unlink(output_file), add = TRUE)

run_multiCSN(
  matrix = matrix(seq_len(8), nrow = 4),
  pseudotime = pseudotime,
  filename = output_file,
  verbose = FALSE
)

stopifnot(identical(captured_pseudotime, pseudotime))
cat("Multi-branch runner regression test passed\n")
