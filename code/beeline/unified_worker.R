#!/usr/bin/env Rscript
# One persistent R process per method; timing and network share the same call.
repo <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(repo, "code", "functions", "functions.R"))
benchmark_root_dir <- local({
  root <- file.path(repo, "code")
  function() root
})
for (method_name in c("GENIE3", "SCORPION", "PIDC", "LEAP", "PPCOR", "SCODE", "SINCERITIES")) {
  source(file.path(repo, "code", "methods", paste0(method_name, ".R")))
}

read_benchmark_case <- function(case_dir) {
  list(
    expression = read.csv(file.path(case_dir, "ExpressionData.csv"),
                          row.names = 1, check.names = FALSE),
    pseudotime = read.csv(file.path(case_dir, "PseudoTime.csv"),
                          row.names = 1, check.names = FALSE)
  )
}

# Run each method using transposed expression matrix convention
method_runners <- list(
  # GENIE3 and PIDC expect genes x cells
  GENIE3 = function(expr, pseudotime, cores, out_file) {
    run_GENIE3(matrix = t(expr), cores = cores, filename = out_file)
  },
  # SCORPION expects genes x cells
  SCORPION = function(expr, pseudotime, cores, out_file) {
    run_SCORPION(matrix = t(expr), n_cores = cores, filename = out_file)
  },
  PIDC = function(expr, pseudotime, cores, out_file) {
    run_PIDC(matrix = t(expr), cores = cores, filename = out_file)
  },
  # LEAP and PPCOR expect genes x cells
  LEAP = function(expr, pseudotime, cores, out_file) {
    run_LEAP(matrix = t(expr), filename = out_file)
  },
  PPCOR = function(expr, pseudotime, cores, out_file) {
    run_PPCOR(matrix = t(expr), filename = out_file)
  },
  # SCODE expects genes x cells
  SCODE = function(expr, pseudotime, cores, out_file) {
    run_SCODE(matrix = t(expr), pseudotime = pseudotime, filename = out_file)
  },
  # SINCERITIES expects genes x cells
  SINCERITIES = function(expr, pseudotime, cores, out_file) {
    run_SINCERITIES(
      matrix = t(expr), pseudotime = pseudotime, cores = cores, filename = out_file
    )
  }
)
source("code/methods/multiCSN.R")

load_average_1min <- function() {
  value <- tryCatch(readLines("/proc/loadavg", n = 1L, warn = FALSE), error = function(e) character())
  if (length(value) != 1L) return(NA_real_)
  suppressWarnings(as.numeric(strsplit(value, " ", fixed = TRUE)[[1L]][1L]))
}

run_benchmark_method <- function(output_root, method, limit = 0L) {
  output_root <- normalizePath(output_root, mustWork = TRUE)
  data.table::setDTthreads(1L)
  runners <- c(list(multiCSN = function(expr, pseudotime, cores, out_file) {
    run_multiCSN(expr, pseudotime, cores = cores, filename = out_file)
  }), method_runners)
  if (method == "multiCSN") loadNamespace("inferCSN")
  paths <- sort(unlist(lapply(c("Curated", "Synthetic"), function(group)
    list.files(file.path(repo, "data/beeline", group), "^ExpressionData\\.csv$",
               recursive = TRUE, full.names = TRUE))))
  if (limit > 0L) paths <- head(paths, limit)
  records <- list()
  for (path in paths) {
    d <- dirname(path)
    group <- basename(dirname(d))
    type <- tolower(basename(dirname(dirname(d))))
    relative <- file.path("networks", paste0(method, "__", basename(d), ".csv"))
    output <- file.path(output_root, relative)
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    load_started <- proc.time()[["elapsed"]]
    bundle <- read_benchmark_case(d)
    expression <- t(as.matrix(bundle$expression))
    load_seconds <- proc.time()[["elapsed"]] - load_started
    set.seed(1L)
    started <- proc.time()[["elapsed"]]
    network <- runners[[method]](expression, bundle$pseudotime, 1L, output)
    seconds <- proc.time()[["elapsed"]] - started
    records[[length(records) + 1L]] <- data.frame(
      Case = basename(d), Method = method, DatasetType = type, DatasetGroup = group,
      RuntimeSeconds = seconds, LoadSeconds = load_seconds, LoadAverage = load_average_1min(),
      Network = relative,
      TimingScope = "inference_and_export", Seed = 1L, Status = "complete")
    data.table::fwrite(do.call(rbind, records),
                      file.path(output_root, "work", paste0(method, ".tsv")), sep = "\t")
    message(length(records), "/", length(paths), " ", method, " ", basename(d), " ", seconds, "s")
  }

}
