# Settings applied by R tasks before loading their analysis packages.
configure_task_runtime <- function(threads = 1L) {
  Sys.setenv(LC_ALL = "C")
  settings <- c(
    OMP_NUM_THREADS = threads, OPENBLAS_NUM_THREADS = threads,
    MKL_NUM_THREADS = threads, VECLIB_MAXIMUM_THREADS = threads,
    NUMEXPR_NUM_THREADS = threads, BLAS_NUM_THREADS = threads,
    BLIS_NUM_THREADS = threads
  )
  do.call(Sys.setenv, as.list(stats::setNames(as.character(settings), names(settings))))
  invisible(NULL)
}

# Build a quoted function invocation for a fresh R process.
r_function_expression <- function(script, name, settings = list()) {
  call <- as.call(c(list(as.name(name)), settings))
  paste0("source(", encodeString(normalizePath(script, mustWork = TRUE), quote = '"'),
         "); ", paste(deparse(call, width.cutoff = 500L), collapse = " "))
}
