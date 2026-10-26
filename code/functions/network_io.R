source("code/functions/table_io.R")

read_fit_checkpoint <- function(path) {
  # gzip is available on the supported macOS/Linux hosts; no R.utils dependency.
  if (grepl("[.]gz$", path)) {
    if (!nzchar(Sys.which("gzip"))) stop("gzip is required for checkpoint reads.")
    return(data.table::fread(cmd = paste("gzip -dc", shQuote(path))))
  }
  data.table::fread(path)
}
