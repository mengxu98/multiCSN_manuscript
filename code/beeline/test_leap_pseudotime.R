suppressPackageStartupMessages(library(data.table))
source("code/methods/LEAP.R")
set.seed(4181)
expression <- matrix(rnorm(6L * 24L), nrow = 6L,
                     dimnames = list(paste0("g", 1:6), paste0("c", sprintf("%02d", 1:24))))
times <- matrix(NA_real_, 24L, 3L, dimnames = list(colnames(expression), c("a", "b", "duplicate_a")))
times[1:16, 1] <- rep(sample(1:8), each = 2L)
times[9:24, 2] <- sample(seq_len(16L))
times[, 3] <- times[, 1]
expected <- rbindlist(lapply(seq_len(ncol(times)), function(branch) {
  cells <- which(is.finite(times[, branch]))
  cells <- cells[order(times[cells, branch], colnames(expression)[cells], method = "radix")]
  answer <- LEAP::MAC_counter(expression[, cells], MAC_cutoff = 0, lag_matrix = FALSE, symmetric = FALSE)
  data.table(regulator = rownames(expression)[answer[, "Row gene index"]],
             target = rownames(expression)[answer[, "Column gene index"]],
             weight = answer[, "Correlation"], branch = branch)
}))
setorderv(expected, c("regulator", "target", "branch"), c(1L, 1L, 1L))
expected <- expected[, .SD[which.max(abs(weight))], by = .(regulator, target)]
expected <- expected[, .(regulator, target, weight)]
canonical <- function(x) {
  z <- as.data.table(x)
  setorder(z, regulator, target)
  z
}
folder <- tempfile("leap-order-check-")
dir.create(folder)
settings <- list(matrix = expression, filename = file.path(folder, "first.csv"))
if ("pseudotime" %in% names(formals(run_LEAP))) settings$pseudotime <- times
observed <- do.call(run_LEAP, settings)
stopifnot("LEAP must agree with branchwise native inference on pseudotime-sorted cells" =
            isTRUE(all.equal(canonical(observed), canonical(expected), tolerance = 1e-12)))
permutation <- sample(seq_len(ncol(expression)))
settings$matrix <- expression[, permutation]
settings$pseudotime <- times[sample(seq_len(nrow(times))), , drop = FALSE]
settings$filename <- file.path(folder, "permuted.csv")
permuted <- do.call(run_LEAP, settings)
stopifnot("Cell storage order must not change pseudotime-ordered LEAP inference" =
            isTRUE(all.equal(canonical(permuted), canonical(expected), tolerance = 1e-12)))
settings$regulators <- c("g1", "g3")
settings$targets <- c("g2", "g5")
settings$filename <- file.path(folder, "filtered.csv")
filtered <- do.call(run_LEAP, settings)
stopifnot(all(filtered$regulator %in% settings$regulators), all(filtered$target %in% settings$targets))
settings$regulators <- NULL
settings$targets <- NULL
settings$pseudotime <- times
rownames(settings$pseudotime)[1] <- "unknown_cell"
stopifnot(inherits(try(do.call(run_LEAP, settings), silent = TRUE), "try-error"))
message("PASS: native branchwise equivalence, tied-state ordering, overlapping and duplicate branches, storage-order invariance, candidate filters and cell identifiers")
