source("code/gse194122/run_network_batch.R")

expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

fixture <- function(count = 4L, fail = FALSE) {
  root <- tempfile("network batch's ")
  dir.create(root)
  source <- file.path(root, "source_manifests")
  dir.create(source)
  tasks <- data.frame(task_order = seq_len(count), task_slug = paste0("task", seq_len(count)),
    state = "MK/E prog", expected_cells = 800L, cell_manifest = paste0("cells", seq_len(count), ".tsv"))
  write.table(tasks, file.path(source, "plan.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  for (name in tasks$cell_manifest) writeLines("cell_id\ncell1", file.path(source, name))
  for (name in c("input.rds", "validation.tsv")) writeLines("fixture", file.path(root, name))
  worker <- file.path(root, "worker.R")
  writeLines(c(
    'fit_state_networks <- function(...) {',
    'args <- list(...)',
    'root <- dirname(dirname(args$output)); slug <- basename(args$output)',
    'cat(slug, Sys.getpid(), as.numeric(Sys.time()), "start", sep="\\t", file=file.path(root,"events.tsv"),append=TRUE)',
    'cat("\\n", file=file.path(root,"events.tsv"),append=TRUE)',
    'stopifnot(Sys.getenv("OMP_NUM_THREADS")=="1", Sys.getenv("LC_ALL")=="C")',
    if (fail) 'if (slug=="task1") {Sys.sleep(0.2); quit(status=7L)} else Sys.sleep(10)' else 'Sys.sleep(0.4)',
    'dir.create(args$output)',
    'saveRDS(args, file.path(args$output,"identity.rds"))',
    'cat(slug, Sys.getpid(), as.numeric(Sys.time()), "end", sep="\\t", file=file.path(root,"events.tsv"),append=TRUE)',
    'cat("\\n", file=file.path(root,"events.tsv"),append=TRUE)',
    '}'
  ), worker)
  list(root = file.path(root, "run"), manifest_source = source,
    input = file.path(root, "input.rds"), validation = file.path(root, "validation.tsv"),
    plan = "plan.tsv", expected_tasks = count, max_jobs = 2L, task_start = 1L, task_end = count,
    response_chunk = 1024L, memory_gate_kib = 1, worker_script = worker)
}

validate_fixture <- function(output, expected) {
  args <- readRDS(file.path(output, "identity.rds"))
  if (!isTRUE(all.equal(args[names(expected)], expected))) stop("Run identity mismatch")
}

run_tests <- function() {
  roots <- character()
  new_fixture <- function(...) {
    settings <- fixture(...)
    roots <<- c(roots, dirname(settings$root))
    settings
  }
  on.exit(unlink(roots, recursive = TRUE), add = TRUE)
  stopifnot(TASKS$donor_deletion$expected_tasks == 49L, TASKS$donor_deletion$max_jobs == 3L,
            TASKS$cell_budget$expected_tasks == 50L, TASKS$cell_budget$max_jobs == 5L,
            TASKS$cell_budget$memory_gate_kib == 125829120, TASKS$cell_budget$verify_exports)
  settings <- new_fixture()
  reads <- 0L
  memory <- function() {reads <<- reads + 1L; if (reads <= 2L) 0 else 2}
  run_batch(settings, memory, validate_fixture)
  stopifnot(file.exists(file.path(settings$root, "run_status/memory_gate.tsv")), reads >= 6L,
            identical(directory_hashes(settings$manifest_source), directory_hashes(file.path(settings$root, "manifests"))))
  events <- read.delim(file.path(settings$root, "events.tsv"), header = FALSE,
                       col.names = c("slug", "pid", "time", "event"))
  events <- events[order(events$time), ]
  active <- cumsum(ifelse(events$event == "start", 1L, -1L))
  stopifnot(nrow(events) == 8L, max(active) == 2L, min(active) >= 0L, tail(active, 1L) == 0L)
  validations <- character()
  check <- function(output, expected) {
    validations <<- c(validations, basename(output)); validate_fixture(output, expected)
  }
  run_batch(modifyList(settings, list(task_end = 1L)), function() stop("Reused outputs must not wait for memory"), check)
  stopifnot(all(paste0("task", 1:4) %in% validations), nrow(read.delim(file.path(settings$root, "events.tsv"), header = FALSE)) == 8L)
  saveRDS(list(), file.path(settings$root, "results/task4/identity.rds"))
  expect_error(run_batch(modifyList(settings, list(validate_only = TRUE)), validate = validate_fixture), "identity mismatch")

  settings <- new_fixture(2L)
  tasks <- read.delim(file.path(settings$manifest_source, "plan.tsv"))
  tasks$task_slug[2L] <- tasks$task_slug[1L]
  write.table(tasks, file.path(settings$manifest_source, "plan.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(run_batch(settings, validate = validate_fixture), "unique directory names")
  stopifnot(!file.exists(file.path(settings$root, "events.tsv")))

  settings <- new_fixture(2L)
  dir.create(settings$root)
  file.copy(settings$manifest_source, settings$root, recursive = TRUE)
  file.rename(file.path(settings$root, "source_manifests"), file.path(settings$root, "manifests"))
  writeLines("changed", file.path(settings$manifest_source, "cells1.tsv"))
  expect_error(run_batch(settings, validate = validate_fixture), "differ from the source")

  settings <- new_fixture(2L)
  file.symlink(settings$input, file.path(settings$manifest_source, "linked.tsv"))
  expect_error(run_batch(settings, validate = validate_fixture), "symbolic link")

  settings <- new_fixture(2L, fail = TRUE)
  start <- Sys.time()
  expect_error(run_batch(settings, function() 2, validate_fixture), "exit 7")
  stopifnot(as.numeric(difftime(Sys.time(), start, units = "secs")) < 5,
    read.delim(file.path(settings$root, "run_status/task1.tsv"))$exit_code == 7L,
    read.delim(file.path(settings$root, "run_status/task2.tsv"))$exit_code == 130L)
  pids <- read.delim(file.path(settings$root, "events.tsv"), header = FALSE)$V2
  stopifnot(!any(pids %in% ps::ps_pids()))
  cat("PASS R batch scheduling, memory gate, task identity, complete resume validation, manifest integrity and failure cancellation\n")
}
run_tests()
