#!/usr/bin/env Rscript
source("code/functions/runtime.R")
configure_task_runtime()

available_memory_kib <- function() {
  line <- grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE)
  if (length(line) != 1L) stop("Linux MemAvailable is required for the memory gate")
  as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1L]][2L])
}

task_arguments <- function(task, settings) list(
  input = settings$input, state = task$state,
  expected_cells = as.integer(task$expected_cells),
  cell_manifest = file.path(settings$manifest_root, task$cell_manifest),
  response_chunk = settings$response_chunk, profile = "gse194122")

run_batch <- function(settings, memory_available = available_memory_kib,
                      validate = NULL) {
  settings$root <- normalizePath(settings$root, mustWork = FALSE)
  root <- settings$root
  source <- normalizePath(settings$manifest_source, mustWork = TRUE)
  settings$manifest_root <- source
  manifests <- source
  plan <- file.path(manifests, settings$plan)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (!all(file.exists(c(settings$input, plan)))) stop("Missing required input")
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
  dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)
  if (is.null(validate)) {
    source("code/functions/network_contract.R", local = TRUE)
    validate <- validate_network_result
  }
  validate_task <- function(i) {
    output <- file.path(root, "results", tasks$task_slug[i])
    if (!dir.exists(output)) stop("Missing result: ", tasks$task_slug[i])
    validate(output)
  }
  active <- list()
  on.exit({
    for (job in active) {
      if (job$process$is_alive()) job$process$kill_tree()
      job$process$wait()
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
        active[[key]] <- NULL
        if (code != 0L) stop("Network worker failed: ", job$slug, " (exit ", code, ")")
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
            message("Waiting for memory: ", task$task_slug, " (", available, " KiB available)")
            last_gate <- as.numeric(Sys.time())
          }
          break
        }
        arguments <- c(list(dataset = paste0("GSE194122_", task$task_slug), output = output), task_arguments(task, settings))
        job <- list(slug = task$task_slug)
        executable <- if (is.null(settings$rscript)) file.path(R.home("bin"), "Rscript") else settings$rscript
        worker <- if (is.null(settings$worker_script)) "code/functions/network_fit.R" else settings$worker_script
        env <- Sys.getenv()
        env["R_PROFILE_USER"] <- "/dev/null"
        if (!is.null(settings$r_library)) env["R_LIBS_USER"] <- settings$r_library
        job$process <- processx::process$new(executable, c("--vanilla", "-e", r_function_expression(worker, "fit_state_networks", arguments)),
          stdout = "", stderr = "",
          env = env, cleanup_tree = TRUE, supervise = TRUE)
        active[[task$task_slug]] <- job
        queue <- queue[-1L]
      }
      if (length(queue) || length(active)) Sys.sleep(0.1)
    }
  }
  for (i in seq_len(nrow(tasks))) validate_task(i)
  cat("All declared network outputs validated.\n")
  invisible(TRUE)
}

TASKS <- list(
  donor_deletion = list(root = "runs/gse194122/donor_deletion_runs",
    manifest_source = "results/gse194122/sampling",
    input = "data/gse194122/gse194122_candidates.rds",
    plan = "donor_deletion_plan.tsv", expected_tasks = 49L, max_jobs = 3L,
    task_start = 1L, task_end = 49L, response_chunk = 1024L, memory_gate_kib = 125829120),
  cell_budget = list(root = "runs/gse194122/cell_budget_runs",
    manifest_source = "results/gse194122/sampling",
    input = "data/gse194122/gse194122_candidates.rds",
    plan = "cell_budget_sampling_plan.tsv", expected_tasks = 50L, max_jobs = 5L,
    task_start = 1L, task_end = 50L, response_chunk = 1024L, memory_gate_kib = 125829120))

run_donor_deletion <- function() run_batch(TASKS$donor_deletion)
run_cell_budget <- function() run_batch(TASKS$cell_budget)
