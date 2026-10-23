write_tsv <- function(x, path, writer = c("fwrite", "table"), na = "NA") {
  writer <- match.arg(writer)
  if (writer == "table" || (!nrow(x) && grepl("[.]gz$", path))) {
    # Keep base R numeric formatting where required. Header-only gzip files
    # also need a connection: some data.table versions omit their gzip trailer.
    connection <- if (grepl("[.]gz$", path)) gzfile(path, "wt") else path
    on.exit(if (inherits(connection, "connection")) close(connection), add = TRUE)
    utils::write.table(x, connection, sep = "\t", quote = FALSE,
                       row.names = FALSE, na = na)
    return(invisible(NULL))
  }
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = na)
}


publish_result_files <- function(stage, output, prefix = "") {
  paths <- list.files(stage, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  paths <- paths[!file.info(file.path(stage, paths))$isdir]
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  for (relative in paths) {
    destination <- file.path(output, paste0(prefix, relative))
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!file.rename(file.path(stage, relative), destination)) {
      stop("Unable to publish result file: ", destination)
    }
  }
  unlink(stage, recursive = TRUE)
  invisible(output)
}
