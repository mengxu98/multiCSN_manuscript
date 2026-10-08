#!/usr/bin/env Rscript

if (!requireNamespace("remotes", quietly = TRUE)) {
  stop("Install remotes before running this script.", call. = FALSE)
}

remotes::install_github(
  c("mengxu98/inferCSN", "mengxu98/multiCSN"),
  dependencies = NA,
  upgrade = "never"
)
