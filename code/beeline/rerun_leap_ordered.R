args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
output <- normalizePath(args[[1L]], mustWork = TRUE)
shard <- as.integer(args[[2L]])
shards <- as.integer(args[[3L]])
stopifnot(shard >= 1L, shard <= shards)
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
source(file.path(output, "code/functions.R"))
source(file.path(output, "code/LEAP.R"))
worker <- file.path(output, paste0("shard_", shard))
if (dir.exists(worker)) stop("Shard output already exists: ", worker)
dir.create(worker)
dir.create(file.path(worker, "networks"))
writeLines(capture.output(sessionInfo()), file.path(worker, "sessionInfo.txt"))
paths <- sort(unlist(lapply(c("Curated", "Synthetic"), function(group) {
  list.files(file.path("data/beeline", group), "^ExpressionData[.]csv$", recursive = TRUE, full.names = TRUE)
})))
stopifnot(length(paths) == 500L, !anyDuplicated(basename(dirname(paths))))
paths <- paths[(seq_along(paths) - 1L) %% shards == shard - 1L]
rows <- list()
for (path in paths) {
  directory <- dirname(path)
  case <- basename(directory)
  expression <- as.matrix(read.csv(path, row.names = 1L, check.names = FALSE))
  pseudotime <- read.csv(file.path(directory, "PseudoTime.csv"), row.names = 1L, check.names = FALSE)
  truth <- as.data.frame(fread(file.path(directory, "refNetwork.csv")))
  destination <- file.path(worker, "networks", paste0("LEAP__", case, ".csv"))
  set.seed(1L)
  started <- Sys.time()
  network <- run_LEAP(expression, pseudotime = pseudotime, filename = destination)
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  saved <- fread(destination)
  stopifnot(nrow(saved) > 0L, all(is.finite(saved$weight)),
            !anyDuplicated(saved[, .(regulator, target)]),
            all(saved$regulator != saved$target), isTRUE(all.equal(as.data.frame(saved), network, check.attributes = FALSE)))
  score <- fast_standard_metric_summary(as.data.frame(saved), truth)
  stopifnot(all(is.finite(unlist(score))))
  rows[[length(rows) + 1L]] <- c(list(Case = case,
    DatasetType = tolower(basename(dirname(dirname(directory)))),
    DatasetGroup = basename(dirname(directory)), Method = "LEAP", RuntimeSeconds = seconds,
    InputTrajectories = ncol(pseudotime), SelectedEdges = nrow(saved), Status = "complete"), score)
  fwrite(rbindlist(rows), file.path(worker, "case_metrics.tsv"), sep = "\t")
  message(length(rows), "/", length(paths), " ", case)
}
writeLines(paste("complete", length(rows)), file.path(worker, "COMPLETE.txt"))
