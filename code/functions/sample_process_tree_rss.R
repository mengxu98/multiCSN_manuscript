#!/usr/bin/env Rscript
# Linux process-tree RSS sum, in KiB, with the same 0.25-second sampling policy.
process_snapshot <- function() {
  paths <- list.dirs("/proc", recursive = FALSE, full.names = TRUE)
  paths <- paths[grepl("^[0-9]+$", basename(paths))]
  rows <- lapply(paths, function(path) tryCatch({
    stat <- readLines(file.path(path, "stat"), warn = FALSE)
    tail <- strsplit(sub("^.*\\) ", "", stat), " ", fixed = TRUE)[[1L]]
    pages <- as.numeric(strsplit(readLines(file.path(path, "statm"), warn = FALSE), " ", fixed = TRUE)[[1L]][2L])
    c(pid = as.numeric(basename(path)), parent = as.numeric(tail[2L]), pages = pages)
  }, error = function(e) NULL))
  if (!length(Filter(Negate(is.null), rows))) return(matrix(numeric(), 0, 3))
  do.call(rbind, rows)
}

sample_rss <- function(root_pid, directory, interval = 0.25) {
  file.create(file.path(directory, "ready"))
  while (!file.exists(file.path(directory, "go"))) Sys.sleep(0.02)
  page_kib <- as.numeric(system2("getconf", "PAGESIZE", stdout = TRUE)) / 1024
  peak <- count <- samples <- 0
  repeat {
    snapshot <- suppressWarnings(process_snapshot())
    if (!nrow(snapshot) || !root_pid %in% snapshot[, "pid"]) break
    selected <- root_pid
    repeat {
      expanded <- union(selected, snapshot[snapshot[, "parent"] %in% selected, "pid"])
      if (length(expanded) == length(selected)) break
      selected <- expanded
    }
    live <- snapshot[snapshot[, "pid"] %in% selected, , drop = FALSE]
    total <- sum(live[, "pages"]) * page_kib
    if (total > peak) { peak <- total; count <- nrow(live) }
    samples <- samples + 1L
    Sys.sleep(interval)
  }
  fields <- c(process_tree_rss_sum_peak_kb = peak, process_tree_processes_at_peak = count,
              process_tree_rss_samples = samples, process_tree_rss_interval_seconds = interval)
  writeLines(paste(names(fields), fields, sep = "\t"), file.path(directory, "summary.tsv"))
}
