#!/usr/bin/env Rscript
source("code/functions/runtime.R")
configure_task_runtime()

available_memory_kib <- function() {
  line <- grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE)
  if (length(line) != 1L) stop("Linux MemAvailable is required for the memory gate")
  as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1L]][2L])
}

directory_hashes <- function(root) {
  paths <- sort(list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                           include.dirs = TRUE, no.. = TRUE), method = "radix")
  if (nzchar(Sys.readlink(root)) || any(nzchar(Sys.readlink(paths))))
    stop("Artifact contains a symbolic link: ", root)
  paths <- paths[!file.info(paths)$isdir]
  stats::setNames(vapply(paths, function(path) digest::digest(file = path, algo = "sha256"), character(1L)),
                  substring(paths, nchar(root) + 2L))
}

task_arguments <- function(task, settings) list(
  input = settings$input, validation = settings$validation, state = task$state,
  expected_cells = as.integer(task$expected_cells),
  cell_manifest = file.path(settings$root, "manifests", task$cell_manifest),
  response_chunk = settings$response_chunk, profile = "gse194122")

verify_network_exports <- function() {
  required <- list(multiCSN = c("NetworkRegions", "Params", "collapse_peak_gene_domains", "factorized_chain_support",
    "find_peaks_near_genes", "get_attribute", "get_peak_annotation", "network_ordinal_weight",
    "parse_peak_ranges", "strongest_mediated_projection"), inferCSN = "fit_greedy_l0_batch")
  for (package in names(required)) {
    missing <- setdiff(required[[package]], getNamespaceExports(package))
    if (length(missing)) stop(package, " missing exports: ", paste(missing, collapse = ", "))
  }
}

run_batch <- function(settings, memory_available = available_memory_kib,
                      validate = NULL) {
  utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  settings$root <- normalizePath(settings$root, mustWork = FALSE)
  root <- settings$root
  source <- normalizePath(settings$manifest_source, mustWork = TRUE)
  manifests <- file.path(root, "manifests")
  plan <- file.path(manifests, settings$plan)
  if (!file.exists(file.path(source, settings$plan))) stop("Missing source task plan")
  hashes <- directory_hashes(source)
  if (dir.exists(manifests)) {
    if (!identical(directory_hashes(manifests), hashes)) stop("Run manifests differ from the source; use a new run directory.")
  } else {
    dir.create(manifests, recursive = TRUE)
    for (relative in names(hashes)) {
      destination <- file.path(manifests, relative)
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(file.path(source, relative), destination)) stop("Manifest copy failed: ", relative)
    }
  }
  if (!identical(directory_hashes(manifests), hashes)) stop("Copied task manifests failed verification")
  # Resolve the newly created tree once it exists, including macOS /var aliases.
  settings$root <- normalizePath(root, mustWork = TRUE)
  root <- settings$root
  if (!all(file.exists(c(settings$input, settings$validation, plan)))) stop("Missing required input")
  tasks <- read.delim(plan, colClasses = "character", check.names = FALSE)
  required <- c("task_order", "task_slug", "state", "expected_cells", "cell_manifest")
  if (!all(required %in% names(tasks)) || nrow(tasks) != settings$expected_tasks || anyNA(tasks[, required]))
    stop("Invalid task plan or task count")
  if (!identical(as.integer(tasks$task_order), seq_len(nrow(tasks)))) stop("Task orders must be contiguous and unique")
  slugs <- tasks$task_slug
  if (anyDuplicated(slugs) || any(!nzchar(slugs) | basename(slugs) != slugs | slugs %in% c(".", "..")))
    stop("Task output names must be unique directory names")
  if (!isTRUE(settings$max_jobs >= 1 && settings$max_jobs <= 5 &&
      settings$task_start >= 1 && settings$task_start <= settings$task_end && settings$task_end <= nrow(tasks)))
    stop("Invalid worker count or task range")
  cells <- tasks$cell_manifest
  if (any(grepl("^/|(^|/)\\.\\.(/|$)", cells)) || !all(file.exists(file.path(manifests, cells))))
    stop("Missing or invalid cell manifest")
  counts <- suppressWarnings(as.integer(tasks$expected_cells))
  if (anyNA(counts) || any(counts < 1L)) stop("Invalid expected cell count")
  for (name in c("results", "logs", "run_status"))
    dir.create(file.path(root, name), recursive = TRUE, showWarnings = FALSE)
  if (is.null(validate)) {
    source("code/functions/network_contract.R", local = TRUE)
    validate <- validate_network_result
  }
  validate_task <- function(i) {
    output <- file.path(root, "results", tasks$task_slug[i])
    if (!dir.exists(output)) stop("Missing result: ", tasks$task_slug[i])
    validate(output, task_arguments(tasks[i, ], settings))
  }
  if (isTRUE(settings$verify_exports)) verify_network_exports()
  active <- list()
  write_status <- function(job, status, code) {
    write.table(data.frame(task = job$slug, status = status, started_utc = job$started,
      finished_utc = if (status == "RUNNING") NA_character_ else utc(), exit_code = code),
      file.path(root, "run_status", paste0(job$slug, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  }
  on.exit({
    for (job in active) {
      if (job$process$is_alive()) job$process$kill_tree()
      job$process$wait()
      write_status(job, "FAIL", 130L)
    }
  }, add = TRUE)
  if (!isTRUE(settings$validate_only)) {
    queue <- seq.int(settings$task_start, settings$task_end)
    last_gate <- -Inf
    while (length(queue) || length(active)) {
      for (key in names(active)) {
        job <- active[[key]]
        if (job$process$is_alive()) next
        code <- job$process$get_exit_status()
        write_status(job, if (code == 0L) "PASS" else "FAIL", code)
        active[[key]] <- NULL
        if (code != 0L) stop("Network worker failed: ", job$slug, " (exit ", code, ")")
        validate_task(job$index)
      }
      while (length(queue) && length(active) < settings$max_jobs) {
        i <- queue[1L]
        task <- tasks[i, ]
        output <- file.path(root, "results", task$task_slug)
        if (file.exists(output)) {
          validate_task(i)
          queue <- queue[-1L]
          next
        }
        available <- memory_available()
        if (!is.finite(available)) stop("Invalid available memory reading")
        if (available < settings$memory_gate_kib) {
          if (as.numeric(Sys.time()) - last_gate >= 30) {
            cat(task$task_slug, available, utc(), sep = "\t", file = file.path(root, "run_status/memory_gate.tsv"), append = TRUE)
            cat("\n", file = file.path(root, "run_status/memory_gate.tsv"), append = TRUE)
            last_gate <- as.numeric(Sys.time())
          }
          break
        }
        arguments <- c(list(dataset = paste0("GSE194122_", task$task_slug), output = output), task_arguments(task, settings))
        job <- list(slug = task$task_slug, index = i, started = utc())
        write_status(job, "RUNNING", NA_integer_)
        executable <- if (is.null(settings$rscript)) file.path(R.home("bin"), "Rscript") else settings$rscript
        worker <- if (is.null(settings$worker_script)) "code/functions/network_fit.R" else settings$worker_script
        env <- Sys.getenv()
        env["R_PROFILE_USER"] <- "/dev/null"
        if (!is.null(settings$r_library)) env["R_LIBS_USER"] <- settings$r_library
        job$process <- tryCatch(processx::process$new(executable, c("--vanilla", "-e", r_function_expression(worker, "fit_state_networks", arguments)),
          stdout = file.path(root, "logs", paste0(task$task_slug, ".log")), stderr = "2>&1",
          env = env, cleanup_tree = TRUE, supervise = TRUE), error = function(e) {
            write_status(job, "FAIL", 130L)
            stop(e)
          })
        active[[task$task_slug]] <- job
        queue <- queue[-1L]
      }
      if (length(queue) || length(active)) Sys.sleep(0.1)
    }
  }
  # A partial rerun must still validate every declared output.
  for (i in seq_len(nrow(tasks))) validate_task(i)
  cat("All declared network outputs validated.\n")
  invisible(TRUE)
}

TASKS <- list(
  donor_deletion = list(root = "runs/gse194122/donor_deletion_runs",
    manifest_source = "results/gse194122/donor_deletion_manifests",
    input = "data/gse194122/gse194122_candidates.rds", validation = "data/gse194122/candidate_validation.tsv",
    plan = "donor_deletion_plan.tsv", expected_tasks = 49L, max_jobs = 3L,
    task_start = 1L, task_end = 49L, response_chunk = 1024L, memory_gate_kib = 125829120),
  cell_budget = list(root = "runs/gse194122/cell_budget_runs",
    manifest_source = "results/gse194122/cell_budget_manifests",
    input = "data/gse194122/gse194122_candidates.rds", validation = "data/gse194122/candidate_validation.tsv",
    plan = "cell_budget_sampling_plan.tsv", expected_tasks = 50L, max_jobs = 5L,
    task_start = 1L, task_end = 50L, response_chunk = 1024L, memory_gate_kib = 125829120,
    verify_exports = TRUE))

run_donor_deletion <- function() run_batch(TASKS$donor_deletion)
run_cell_budget <- function() run_batch(TASKS$cell_budget)
