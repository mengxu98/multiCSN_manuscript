#!/usr/bin/env Rscript

if (!requireNamespace("remotes", quietly = TRUE)) {
  stop("Install remotes before running this script.", call. = FALSE)
}

refs <- c(
  inferCSN = "mengxu98/inferCSN@03010a34e6b7247e34b52d1feb83be007cd10cd3",
  multiCSN = "mengxu98/multiCSN@984d0fe75e2ead58aec6955ea60ac84aa6ddaf5e"
)
for (ref in refs) {
  remotes::install_github(ref, dependencies = NA, upgrade = "never")
}
for (pkg in names(refs)) {
  expected <- sub(".*@", "", refs[[pkg]])
  observed <- utils::packageDescription(pkg)[["RemoteSha"]]
  if (!identical(observed, expected)) stop("Commit verification failed for ", pkg)
}
