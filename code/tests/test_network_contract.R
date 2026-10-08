suppressPackageStartupMessages(library(data.table))
source("code/functions/network_io.R")
source("code/functions/network_contract.R")
root <- tempfile("network-contract-")
dir.create(root)
files <- c("tf_gene.tsv", "tf_region.tsv", "region_gene.tsv", "triplets.tsv",
           "mediated_tf_gene.tsv", "tf_gene_accounting.tsv", "tf_region_accounting.tsv",
           "region_gene_accounting.tsv")
for (file in files) fwrite(data.table(regulator = "TF", target = "G", region = "chr1:1-2"),
                          file.path(root, file), sep = "\t")
fwrite(data.table(endpoint = c("TF-gene", "TF-region", "region-gene", "factorized-chain",
                              "mediated-TF-gene"), selected = 1L),
       file.path(root, "endpoint_summary.tsv"), sep = "\t")
fails <- function(expr) inherits(tryCatch({force(expr); NULL}, error = base::identity), "error")
stopifnot(validate_network_result(root), fails(validate_network_result(root, list(state = "HSC"))))
input_root <- tempfile("network-contract-inputs-")
dir.create(input_root)
input_files <- file.path(input_root, c("input.rds", "validation.tsv", "cells.tsv"))
for (path in input_files) writeLines("fixture", path)
input_files <- normalizePath(input_files, mustWork = TRUE)
inputs <- normalizePath(c("code/functions/network_contract.R", input_files))
outputs <- list.files(root, full.names = TRUE)
settings <- list(input = input_files[1L], validation = input_files[2L],
  state = "HSC", expected_cells = 10L, cell_manifest = input_files[3L],
  response_chunk = 1024L, profile = "gse194122")
identity <- list(settings = settings,
  inputs = data.table(file = inputs, sha256 = unname(sha256_files(inputs))),
  outputs = data.table(file = basename(outputs), sha256 = unname(sha256_files(outputs))))
saveRDS(identity, file.path(root, "run_provenance.rds"))
stopifnot(validate_network_result(root, list(state = "HSC", expected_cells = 10L)),
          fails(validate_network_result(root, list(state = "Normoblast"))),
          fails(validate_network_result(root, list(expected_cells = 9L))))
run_validator <- function(settings) {
  invocation <- as.call(c(list(as.name("validate_network_result")), list(root = root, expected = settings)))
  expression <- paste('source("code/functions/network_contract.R")',
    paste(deparse(invocation, width.cutoff = 500L), collapse = " "),
    'cat("PASS verified network inputs")', sep = "; ")
  output <- suppressWarnings(system2(Sys.which("Rscript"),
    c("--vanilla", "-e", shQuote(expression)), stdout = TRUE, stderr = TRUE))
  list(status = if (is.null(attr(output, "status"))) 0L else attr(output, "status"),
       output = paste(output, collapse = "\n"))
}
valid <- run_validator(settings)
wrong_state <- run_validator(modifyList(settings, list(state = "Normoblast")))
stopifnot(valid$status == 0L,
          grepl("PASS verified network inputs", valid$output, fixed = TRUE),
          wrong_state$status != 0L,
          grepl("Run identity mismatch: state", wrong_state$output, fixed = TRUE))
fwrite(data.table(regulator = "OTHER", target = "G", region = "chr1:1-2"),
       file.path(root, "tf_gene.tsv"), sep = "\t")
stopifnot(fails(validate_network_result(root, list(state = "HSC"))))
corrupt <- run_validator(settings)
stopifnot(corrupt$status != 0L,
          grepl("Output content changed", corrupt$output, fixed = TRUE))
unlink(file.path(root, "triplets.tsv"))
stopifnot(fails(validate_network_result(root)))
unlink(root, recursive = TRUE)
unlink(input_root, recursive = TRUE)
cat("PASS shared function and subprocess validation: identity, state/cell count, content and completeness checked\n")
