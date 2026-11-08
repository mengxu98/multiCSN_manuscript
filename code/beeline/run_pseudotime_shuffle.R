args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
output <- args[[1L]]
if (file.exists(output)) stop("Output already exists: ", output)
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
source("code/functions/functions.R")
source("code/methods/multiCSN.R")
dir.create(output, recursive = TRUE)
dir.create(file.path(output, "networks"))
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
paths <- sort(unlist(lapply(c("Curated", "Synthetic"), function(group) {
  list.files(file.path("data/beeline", group), "^ExpressionData[.]csv$", recursive = TRUE, full.names = TRUE)
})))
stopifnot(length(paths) == 500L, !anyDuplicated(basename(dirname(paths))))
seeds <- c(61105L, 61106L, 61107L)
rows <- list()
for (index in seq_along(paths)) {
  path <- paths[[index]]
  directory <- dirname(path)
  case <- basename(directory)
  expression <- t(as.matrix(read.csv(path, row.names = 1L, check.names = FALSE)))
  times <- read.csv(file.path(directory, "PseudoTime.csv"), row.names = 1L, check.names = FALSE)
  truth <- as.data.frame(fread(file.path(directory, "refNetwork.csv")))
  for (seed in seeds) {
    set.seed(seed)
    shuffled <- times
    for (branch in seq_len(ncol(times))) {
      previous <- which(vapply(seq_len(branch - 1L), function(j) identical(times[[j]], times[[branch]]), logical(1)))
      if (length(previous)) {
        shuffled[[branch]] <- shuffled[[previous[[1L]]]]
      } else {
        states <- sort(unique(times[[branch]][is.finite(times[[branch]])]))
        replacement <- states[sample.int(length(states))]
        finite <- is.finite(times[[branch]])
        shuffled[finite, branch] <- replacement[match(times[finite, branch], states)]
      }
      stopifnot(identical(is.na(shuffled[[branch]]), is.na(times[[branch]])))
      finite <- is.finite(times[[branch]])
      original_groups <- match(times[finite, branch], unique(times[finite, branch]))
      shuffled_groups <- match(shuffled[finite, branch], unique(shuffled[finite, branch]))
      stopifnot(identical(original_groups, shuffled_groups))
    }
    destination <- file.path(output, "networks", paste0(case, "__", seed, ".csv"))
    network <- run_multiCSN(expression, shuffled, cores = 1L, filename = destination)
    saved <- fread(destination)
    stopifnot(is.data.frame(network), all(is.finite(saved$weight)),
              !anyDuplicated(saved[, .(regulator, target)]))
    score <- fast_standard_metric_summary(as.data.frame(saved), truth)
    stopifnot(all(is.finite(unlist(score))))
    rows[[length(rows) + 1L]] <- c(list(Case = case,
      DatasetType = tolower(basename(dirname(dirname(directory)))),
      DatasetGroup = basename(dirname(directory)), ShuffleSeed = seed,
      OrderChanged = !identical(shuffled, times), Status = "complete"), score)
  }
  if (index %% 25L == 0L) {
    fwrite(rbindlist(rows), file.path(output, "case_metrics.tsv"), sep = "\t")
    message(index, "/500 cases; ", length(rows), " shuffled fits")
  }
}
stopifnot(length(rows) == 1500L)
fwrite(rbindlist(rows), file.path(output, "case_metrics.tsv"), sep = "\t")
writeLines("complete: 500 cases x 3 fixed shuffle seeds; original state memberships and branch masks preserved", file.path(output, "COMPLETE.txt"))
