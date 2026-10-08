suppressPackageStartupMessages(library(data.table))
source("code/functions/table_io.R")

read_function <- function(path, name, environment) {
  expressions <- parse(path)
  for (expression in expressions) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        identical(expression[[2L]], as.name(name))) {
      eval(expression, environment)
      return(environment[[name]])
    }
  }
  stop("Missing function: ", name)
}

root <- tempfile("linger-endpoints-")
dir.create(root)
network_root <- file.path(root, "raw")
endpoint_root <- file.path(root, "endpoints")
dir.create(file.path(endpoint_root, "Pando"), recursive = TRUE)
writeLines("keep", file.path(endpoint_root, "Pando", "existing.txt"))
timepoints <- c("day 7", "day 9", "day 11", "day 14")
targets <- paste0("TF", 1:19)
for (timepoint in timepoints) {
  folder <- file.path(network_root, gsub(" ", "_", timepoint), "networks")
  dir.create(folder, recursive = TRUE)
  write_tsv(data.table(region = "chr1-1-10", regulator = targets, score = 1),
            file.path(folder, "tf_region.tsv.gz"))
  write_tsv(data.table(region = "chr1-1-10", target = "G1", score = 1),
            file.path(folder, "region_gene.tsv.gz"))
  write_tsv(data.table(endpoint = "region_gene", relations = 1),
            file.path(folder, "network_summary.tsv"))
}
path <- "code/gse274113/evaluate_linger_multilayer.R"
prepare <- read_function(path, "prepare_endpoints", environment())
prepare()
prepare()
stopifnot(readLines(file.path(endpoint_root, "Pando", "existing.txt")) == "keep",
          nrow(fread(file.path(endpoint_root, "LINGER", "endpoint_summary.tsv"))) == 4L)

# Dispatch passes the selected paths and method into the evaluator's real
# function interface for isolated inputs and outputs.
dispatch_environment <- new.env(parent = globalenv())
dispatch_environment$arguments <- list(evidence = root, object = root)
dispatch_environment$endpoint_root <- endpoint_root
dispatch_environment$truth_root <- root
dispatch_environment$output_root <- file.path(root, "evaluation")
captured <- NULL
dispatch_environment$sys.source <- function(file, envir) {
  name <- sub("[.]R$", "", basename(file))
  envir[[name]] <- function(arguments) captured <<- list(name = name, arguments = arguments)
}
dispatch <- read_function(path, "run_endpoint_evaluator_for_linger", dispatch_environment)
for (mode in c("binding", "contact")) {
  dispatch(mode)
  stopifnot(captured$arguments$methods == "LINGER",
            captured$arguments$endpoints == endpoint_root,
            captured$arguments$output == dispatch_environment$output_root,
            captured$name == if (mode == "binding") "evaluate_tf_region_binding" else
              "evaluate_region_gene_contact")
  evaluator <- read_function(file.path("code/gse274113", paste0(captured$name, ".R")),
                             captured$name, new.env(parent = globalenv()))
  stopifnot("arguments" %in% names(formals(evaluator)))
}
unlink(root, recursive = TRUE)
cat("PASS LINGER endpoint publication and binding/contact function dispatch\n")
