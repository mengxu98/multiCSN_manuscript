source("code/functions/table_io.R")

sha256_files <- function(paths, batch_size = 256L) {
  paths <- normalizePath(paths, mustWork = TRUE)
  if (!length(paths)) {
    return(character())
  }
  if (anyDuplicated(paths)) stop("SHA-256 input paths must be unique.")
  executable <- if (nzchar(Sys.which("sha256sum"))) "sha256sum" else "shasum"
  starts <- seq.int(1L, length(paths), by = as.integer(batch_size))
  result <- character(length(paths))
  for (first in starts) {
    indices <- first:min(length(paths), first + batch_size - 1L)
    arguments <- if (executable == "shasum") {
      c("-a", "256", shQuote(paths[indices]))
    } else {
      shQuote(paths[indices])
    }
    output <- system2(executable, arguments, stdout = TRUE, stderr = TRUE)
    if (length(output) != length(indices) || any(!grepl("^[0-9a-f]{64}[[:space:]]", output))) {
      stop("Batched SHA-256 failed for paths ", first, " through ", max(indices), ".")
    }
    result[indices] <- substr(output, 1L, 64L)
  }
  names(result) <- paths
  result
}

read_fit_checkpoint <- function(path) {
  # gzip is available on the supported macOS/Linux hosts; no R.utils dependency.
  if (grepl("[.]gz$", path)) {
    if (!nzchar(Sys.which("gzip"))) stop("gzip is required for checkpoint reads.")
    return(data.table::fread(cmd = paste("gzip -dc", shQuote(path))))
  }
  data.table::fread(path)
}
