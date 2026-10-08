source("code/functions/runtime.R")

#!/usr/bin/env Rscript

# R is the only scaling entry. Each measured case runs in a fresh R process.
prepare_scaling_inputs <- function(repo, output_root) {
  manifest_path <- file.path(repo, "data/beeline/scaling/scaling_manifest.csv")
  manifest <- read.csv(manifest_path)
  dir.create(file.path(output_root, "inputs"), recursive = TRUE, showWarnings = FALSE)
  stopifnot(nrow(manifest) == 10L)
  for (i in seq_len(nrow(manifest))) {
    m <- manifest[i, , drop = FALSE]
    path <- file.path(output_root, "inputs", paste0(m$Case, ".rds"))
    d <- file.path(repo, "data/beeline/scaling", m$Scenario, m$Case)
    source_files <- c(manifest_path, file.path(d, c("ExpressionData.csv", "PseudoTime.csv")))
    if (!all(file.exists(source_files))) stop("Missing scaling input: ", m$Case)
    x <- t(as.matrix(read.csv(file.path(d, "ExpressionData.csv"), row.names = 1L,
                              check.names = FALSE)))
    pt <- read.csv(file.path(d, "PseudoTime.csv"), row.names = 1L,
                   check.names = FALSE)[[1L]]
    stopifnot(nrow(x) == m$RequestedCells, ncol(x) == m$TotalGeneCount,
              length(pt) == nrow(x), all(is.finite(pt)))
    saveRDS(list(x = x, pt = pt, metadata = m), path, compress = FALSE)
  }
}

measure_scaling_case <- function(input_file, output_dir, cores) {
  suppressPackageStartupMessages(library(inferCSN))
  z <- readRDS(input_file)
  gc()
  file.create(file.path(output_dir, "ready"))
  while (!file.exists(file.path(output_dir, "go"))) Sys.sleep(0.01)
  clock <- proc.time()
  network <- inferCSN::inferCSN(z$x, pseudotime = z$pt, cores = cores, verbose = FALSE)
  timing <- proc.time() - clock
  file.create(file.path(output_dir, "done"))
  saveRDS(network, file.path(output_dir, "network.rds"), compress = FALSE)
  m <- z$metadata
  result <- data.frame(
    Case = m$Case, Scenario = m$Scenario,
    AxisValue = if (m$Scenario == "cell_scaling") nrow(z$x) else ncol(z$x),
    AxisLabel = if (m$Scenario == "cell_scaling") "Cells" else "Total genes",
    RequestedGenes = m$RequestedGenes, RequestedTFs = m$RequestedTFs,
    Cores = cores, RequestedCores = cores,
    RuntimeSeconds = unname(timing[["elapsed"]]),
    UserSeconds = unname(timing[["user.self"]] + timing[["user.child"]]),
    SystemSeconds = unname(timing[["sys.self"]] + timing[["sys.child"]]),
    InputCells = nrow(z$x), InputGenes = ncol(z$x), OutputEdges = nrow(network),
    TimingScope = "public_inference_only"
  )
  write.csv(result, file.path(output_dir, "timing.csv"), row.names = FALSE)
}

scaling_tree_rss <- function(handle) {
  children <- tryCatch(ps::ps_children(handle, recursive = TRUE), error = function(e) {
    if (!ps::ps_is_running(handle)) return(list())
    stop(e)
  })
  handles <- c(list(handle), children)
  rss <- vapply(handles, function(p) {
    tryCatch(unname(ps::ps_memory_info(p)[["rss"]]), error = function(e) {
      if (inherits(e, c("no_such_process", "zombie_process")) || !ps::ps_is_running(p) ||
          identical(tryCatch(ps::ps_status(p), error = function(e) NULL), "zombie")) {
        return(NA_real_)
      }
      stop(e)
    })
  }, numeric(1L))
  list(mib = sum(rss, na.rm = TRUE) / 1024^2, processes = sum(!is.na(rss)))
}

run_scaling_worker <- function(input_file, output_dir, cores, script, rscript, interval) {
  log <- file.path(output_dir, "worker.log")
  process <- processx::process$new(rscript,
    c("--vanilla", "-e", r_function_expression(script, "measure_scaling_case",
      list(input_file = input_file, output_dir = output_dir, cores = as.integer(cores)))),
    stdout = log, stderr = "2>&1", cleanup_tree = TRUE, supervise = TRUE)
  on.exit({
    if (process$is_alive()) process$kill_tree()
    process$wait()
  }, add = TRUE)
  ensure_alive <- function() {
    if (!process$is_alive() && !file.exists(file.path(output_dir, "done"))) {
      stop("Scaling worker failed: ", output_dir, "\n",
           paste(tail(readLines(log, warn = FALSE), 40L), collapse = "\n"))
    }
  }
  while (!file.exists(file.path(output_dir, "ready"))) {
    ensure_alive()
    Sys.sleep(0.05)
  }
  handle <- ps::ps_handle(process$get_pid())
  peak <- baseline <- scaling_tree_rss(handle)
  samples <- 1L
  previous <- unname(proc.time()[["elapsed"]])
  max_interval <- 0
  file.create(file.path(output_dir, "go"))
  while (!file.exists(file.path(output_dir, "done"))) {
    ensure_alive()
    current <- tryCatch(scaling_tree_rss(handle), error = function(e) {
      # A completed worker can exit between the marker check and RSS query.
      if (file.exists(file.path(output_dir, "done"))) return(NULL)
      stop(e)
    })
    if (is.null(current)) break
    if (current$mib > peak$mib) peak <- current
    now <- unname(proc.time()[["elapsed"]])
    max_interval <- max(max_interval, now - previous)
    previous <- now
    samples <- samples + 1L
    Sys.sleep(interval)
  }
  process$wait()
  if (process$get_exit_status() != 0L || !file.exists(file.path(output_dir, "network.rds"))) {
    stop("Scaling output failed: ", output_dir, "\n",
         paste(tail(readLines(log, warn = FALSE), 40L), collapse = "\n"))
  }
  row <- read.csv(file.path(output_dir, "timing.csv"))
  stopifnot(nrow(row) == 1L)
  row$MemoryTreePeakMb <- peak$mib
  row$MemoryTreeBeforeMb <- baseline$mib
  row$MemoryTreeProcessesAtPeak <- peak$processes
  row$MemorySamples <- samples
  row$MemorySamplingIntervalSeconds <- interval
  row$MemoryMaxObservedIntervalSeconds <- max_interval
  row
}

scaling_directory_hashes <- function(root) {
  entries <- list.files(root, recursive = TRUE, all.files = TRUE,
                        include.dirs = TRUE, no.. = TRUE)
  paths <- file.path(root, entries)
  if (any(nzchar(Sys.readlink(c(root, paths))))) stop("Scaling output contains symbolic links")
  files <- sort(entries[!file.info(paths)$isdir])
  stats::setNames(lapply(file.path(root, files), digest::digest,
                         algo = "sha256", file = TRUE), files)
}

publish_scaling_result <- function(candidate, target, backup_root) {
  if (nzchar(Sys.readlink(target))) stop("Scaling output must be a real directory")
  if (!dir.exists(target)) {
    if (!file.rename(candidate, target)) stop("Cannot publish scaling result")
    return(invisible(NULL))
  }
  original <- scaling_directory_hashes(target)
  dir.create(backup_root, recursive = TRUE, showWarnings = FALSE)
  backup <- tempfile("previous-", tmpdir = backup_root)
  previous <- file.path(backup, "contents")
  dir.create(previous, recursive = TRUE)
  entries <- list.files(target, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (!all(file.copy(entries, previous, recursive = TRUE, copy.date = TRUE)) ||
      !identical(scaling_directory_hashes(previous), original) ||
      !identical(scaling_directory_hashes(target), original)) stop("Scaling backup verification failed")
  jsonlite::write_json(original, file.path(backup, "sha256.json"), auto_unbox = TRUE, pretty = TRUE)
  rollback <- tempfile(".core_scaling-rollback-", tmpdir = dirname(target))
  if (!file.rename(target, rollback)) stop("Cannot move previous scaling result")
  published <- FALSE
  on.exit({
    if (!published && !file.rename(rollback, target)) warning("Restore previous result from: ", rollback)
  }, add = TRUE)
  if (!file.rename(candidate, target)) stop("Cannot publish scaling result")
  published <- TRUE
  unlink(rollback, recursive = TRUE)
  invisible(backup)
}

run_core_scaling <- function(repo = getwd(), output_root = file.path(repo, "results/beeline/core_scaling"),
                             backup_root = path.expand("~/backups/multiCSN_manuscript/core_scaling"),
                             rscript = file.path(R.home("bin"), "Rscript"), interval = 0.05) {
  stopifnot(length(interval) == 1L, is.finite(interval), interval > 0)
  for (package in c("processx", "ps", "digest", "jsonlite")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Missing scaling package: ", package)
  }
  repo <- normalizePath(repo, mustWork = TRUE)
  source(file.path(repo, "code/functions/runtime.R"), local = TRUE)
  configure_task_runtime()
  source(file.path(repo, "code/beeline/summarize_scaling.R"), local = TRUE)
  script <- file.path(repo, "code/beeline/run_core_scaling.R")
  dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
  work <- tempfile(".scaling-", tmpdir = dirname(output_root))
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  prepare_scaling_inputs(repo, work)
  inputs <- sort(list.files(file.path(work, "inputs"), pattern = "\\.rds$", full.names = TRUE))
  stopifnot(length(inputs) == 10L)
  rows <- list()
  for (input in inputs) {
    case <- tools::file_path_sans_ext(basename(input))
    for (cores in c(1L, 2L, 4L, 6L, 8L, 10L)) {
      directory <- file.path(work, "runs", case, paste0("cores_", cores))
      dir.create(directory, recursive = TRUE)
      row <- run_scaling_worker(input, directory, cores, script, rscript, interval)
      rows[[length(rows) + 1L]] <- row
      message(sprintf("SCALING_DONE %d/60 %s cores=%d seconds=%.3f peak_MiB=%.1f",
                      length(rows), case, cores, row$RuntimeSeconds, row$MemoryTreePeakMb))
    }
  }
  write.csv(do.call(rbind, rows), file.path(work, "core_scaling_runs.csv"), row.names = FALSE)
  candidate <- file.path(work, "result")
  summarize_scaling(work, candidate)
  stopifnot(all(read.csv(file.path(candidate, "edge_count_invariance.csv"))$FullOutputExact))
  if (!file.rename(file.path(work, "core_scaling_runs.csv"), file.path(candidate, "core_scaling_runs.csv"))) {
    stop("Cannot finalize scaling run table")
  }
  publish_scaling_result(candidate, output_root, backup_root)
  message("SCALING_COMPLETE 60")
  invisible(output_root)
}

if (sys.nframe() == 0L) {
  configure_task_runtime()
  run_core_scaling()
}
