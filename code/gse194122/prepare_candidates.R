source("code/functions/candidate_preparation.R")

prepare_candidates <- function(input = "data/gse194122/processed_multiome.h5ad",
    destination = "data/gse194122/gse194122_candidates.rds") {
  prepare_dataset_candidates("GSE194122", input, destination)
}

if (sys.nframe() == 0L) prepare_candidates()
