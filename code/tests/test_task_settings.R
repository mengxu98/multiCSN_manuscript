source("code/functions/runtime.R")
check_task_library <- function() {
  previous_wd <- getwd()
  previous_lib <- .libPaths()
  root <- tempfile("task-library-")
  dir.create(file.path(root, "analysis"), recursive = TRUE)
  library <- file.path(root, "comparator_r_library")
  dir.create(library)
  on.exit({
    setwd(previous_wd)
    .libPaths(previous_lib)
    unlink(root, recursive = TRUE)
  }, add = TRUE)
  setwd(file.path(root, "analysis"))
  configure_task_runtime(comparator_library = TRUE)
  stopifnot(identical(.libPaths()[[1L]], normalizePath(library)),
            identical(Sys.getenv("OMP_NUM_THREADS"), "1"))
}
check_task_library()
check_control_output_collisions <- function() {
  root <- tempfile("control-collisions-")
  network_root <- file.path(root, "results/gse274113/networks")
  dir.create(network_root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  guards <- Filter(function(expr) is.call(expr) && identical(expr[[1L]], as.name("if")) &&
    identical(expr[[2L]], quote(is.null(script_options$output))),
    as.list(parse("code/gse274113/build_control_object.R")))
  stopifnot(length(guards) == 1L)
  env <- new.env(parent = baseenv())
  env$repo_root <- root
  env$script_options <- list(output = NULL)
  eval(guards[[1L]], env)
  for (name in c("multiCSN", "Pando", "SCENICplus", "endpoint_summary.tsv", "endpoint_availability.tsv")) {
    path <- file.path(network_root, name)
    writeLines("existing formal result", path)
    error <- tryCatch({eval(guards[[1L]], env); NULL}, error = identity)
    stopifnot(inherits(error, "error"), grepl(path, conditionMessage(error), fixed = TRUE),
              identical(readLines(path), "existing formal result"))
    env$script_options <- list(output = file.path(root, "isolated-control-object"))
    eval(guards[[1L]], env)
    env$script_options <- list(output = NULL)
    unlink(path)
  }
}
check_control_output_collisions()
# Exercise the canonical five-state entry without fitting or writing formal results.
expressions <- as.list(parse("code/gse194122/fit_state_networks.R"))
entry <- Filter(function(expr) is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
  identical(expr[[2L]], as.name("fit_canonical_states")), expressions)
stopifnot(length(entry) == 1L)
env <- new.env()
env$calls <- list()
env$fit_state_networks <- function(...) {
  settings <- list(...)
  env$calls[[length(env$calls) + 1L]] <- settings
  data.table::data.table(cell_id = paste0(settings$state, seq_len(settings$expected_cells)),
                         state = settings$state)
}
env$file.exists <- function(...) FALSE
env$dir.create <- function(...) invisible(TRUE)
env$write_tsv <- function(table, path) {
  stopifnot(nrow(table) == 10952L, !anyDuplicated(table$cell_id))
}
env$rbindlist <- data.table::rbindlist
env$setorder <- data.table::setorder
eval(entry[[1L]], env)
env$fit_canonical_states()
stopifnot(length(env$calls) == 5L,
  identical(vapply(env$calls, `[[`, character(1L), "state"),
            c("HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast")))
cat("PASS task runtime, output collision guards, and five canonical state settings\n")
