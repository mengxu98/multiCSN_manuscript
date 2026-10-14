#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L) {
  snapshot <- normalizePath(args[[1L]], mustWork = TRUE)
  library <- normalizePath(args[[2L]], mustWork = FALSE)
  manifest <- jsonlite::fromJSON(file.path(snapshot, "package_sources.json"))
  selected <- match(c("inferCSN", "multiCSN"), manifest$package)
  if (anyNA(selected) || anyDuplicated(manifest$package)) {
    stop("The snapshot must identify inferCSN and multiCSN exactly once.", call. = FALSE)
  }
  manifest <- manifest[selected, , drop = FALSE]
  archives <- file.path(snapshot, basename(manifest$source_archive))
  hashes <- vapply(archives, digest::digest, character(1), file = TRUE, algo = "sha256")
  if (!identical(unname(hashes), manifest$sha256)) {
    stop("Snapshot source archive SHA-256 mismatch.", call. = FALSE)
  }
  dir.create(library, recursive = TRUE, showWarnings = FALSE)
  library <- normalizePath(library, mustWork = TRUE)
  .libPaths(c(library, .libPaths()))
  for (archive in archives) {
    status <- system2(
      file.path(R.home("bin"), "R"),
      args = c("CMD", "INSTALL", paste0("--library=", shQuote(library)), shQuote(archive)),
      env = paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep)))
    )
    if (status != 0L) {
      stop("Snapshot installation failed: ", basename(archive), call. = FALSE)
    }
  }
} else if (length(args) == 0L) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("Install remotes before running this script.", call. = FALSE)
  }
  remotes::install_github(
    c("mengxu98/inferCSN", "mengxu98/multiCSN"),
    dependencies = NA,
    upgrade = "never"
  )
} else {
  stop("Usage: install_csn_packages.R [snapshot_directory isolated_library]", call. = FALSE)
}
