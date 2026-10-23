source("code/functions/candidate_preparation.R")

prepare_candidates <- function(input = "data/gse162170/release",
    destination = "data/gse162170/gse162170_candidates.rds") {
  prepare_dataset_candidates("GSE162170", input, destination)
}

if (sys.nframe() == 0L) prepare_candidates()
