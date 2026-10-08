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
