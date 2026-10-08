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

sample_rss <- function(root_pid, directory, producer, stage, command, interval = 0.25) {
  utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  hash <- function(path) digest::digest(file = path, algo = "sha256")
  host <- system2("hostname", "-f", stdout = TRUE)
  cpus <- parallel::detectCores()
  model <- grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[1L]
  model <- trimws(sub("^[^:]*:", "", model))
  kernel <- unname(Sys.info()[["release"]])
  machine_hash <- hash("/etc/machine-id")
  fingerprint <- digest::digest(paste0(host, "\t", machine_hash, "\t", cpus, "\t", model, "\t", kernel, "\n"),
                                algo = "sha256", serialize = FALSE)
  run_id <- paste(stage, format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), Sys.getpid(), sep = "_")
  thread_names <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                    "NUMEXPR_NUM_THREADS", "BLAS_NUM_THREADS", "BLIS_NUM_THREADS")
  fields <- c(resource_block_started = run_id, resource_block_started_utc = utc(),
    resource_stage = stage, resource_requested_workers = 1,
    resource_host_name = host, resource_machine_id_sha256 = machine_hash,
    resource_boot_id = readLines("/proc/sys/kernel/random/boot_id"), resource_host_fingerprint = fingerprint,
    resource_cpu_online = cpus, resource_cpu_model = model, resource_kernel_release = kernel,
    resource_producer_path = normalizePath(producer), resource_producer_sha256 = hash(producer),
    resource_wrapper_sha256 = hash("code/functions/runtime.sh"),
    resource_sampler_sha256 = hash("code/functions/sample_process_tree_rss.R"),
    resource_thread_environment = paste(paste(thread_names, Sys.getenv(thread_names), sep = "="), collapse = ";"),
    resource_command = command,
    resource_command_started_epoch_ns = sprintf("%.0f", as.numeric(Sys.time()) * 1e9))
  writeLines(paste(names(fields), fields, sep = "\t"), file.path(directory, "header.tsv"))
  writeLines(run_id, file.path(directory, "run_id"))
  file.create(file.path(directory, "ready"))
  while (!file.exists(file.path(directory, "go"))) Sys.sleep(0.02)
  page_kib <- as.numeric(system2("getconf", "PAGESIZE", stdout = TRUE)) / 1024
  stopifnot(is.finite(page_kib), page_kib > 0)
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
              process_tree_rss_samples = samples, process_tree_rss_interval_seconds = interval,
              resource_command_completed_epoch_ns = sprintf("%.0f", as.numeric(Sys.time()) * 1e9))
  writeLines(paste(names(fields), fields, sep = "\t"), file.path(directory, "summary.tsv"))
}
