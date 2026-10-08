suppressPackageStartupMessages(library(data.table))
source("code/functions/network_io.R")

root <- tempfile("network-io-")
dir.create(root)
table_in <- data.table(target = c("g1", "g2"), value = c(0.5, 1.5))
write_tsv(table_in, file.path(root, "plain.tsv"))
write_tsv(table_in, file.path(root, "compressed.tsv.gz"))
stopifnot(identical(as.data.frame(read_fit_checkpoint(file.path(root, "plain.tsv"))),
                    as.data.frame(table_in)))
stopifnot(identical(as.data.frame(read_fit_checkpoint(file.path(root, "compressed.tsv.gz"))),
                    as.data.frame(table_in)))
hashes <- sha256_files(c(file.path(root, "plain.tsv"), file.path(root, "compressed.tsv.gz")))
stopifnot(length(hashes) == 2L, all(grepl("^[0-9a-f]{64}$", hashes)),
          length(unique(hashes)) == 2L)

# Both established serializers keep their numeric and missing-value formatting.
format_table <- data.frame(
  label = c('quoted "gene"', "missing"),
  value = c(pi, NA_real_), flag = c(TRUE, FALSE)
)
for (writer in c("fwrite", "table")) {
  for (na_value in c("NA", "")) {
    reference <- file.path(root, "reference.tsv")
    actual <- file.path(root, "actual.tsv")
    if (writer == "fwrite") {
      data.table::fwrite(format_table, reference, sep = "\t", quote = FALSE, na = na_value)
    } else {
      utils::write.table(format_table, reference, sep = "\t", quote = FALSE,
                         row.names = FALSE, na = na_value)
    }
    write_tsv(format_table, actual, writer = writer, na = na_value)
    stopifnot(identical(readLines(reference), readLines(actual)))
    compressed <- file.path(root, "format.tsv.gz")
    write_tsv(format_table, compressed, writer = writer, na = na_value)
    stopifnot(identical(readLines(compressed), readLines(reference)))
  }
}

# Empty compressed exports must have a valid gzip trailer and retain headers.
empty <- table_in[0]
empty_path <- file.path(root, "empty.tsv.gz")
write_tsv(empty, empty_path)
stopifnot(identical(names(read_fit_checkpoint(empty_path)), names(empty)),
          nrow(read_fit_checkpoint(empty_path)) == 0L,
          identical(system2("gzip", c("-t", shQuote(empty_path))), 0L))

invalid <- tryCatch(write_tsv(table_in, file.path(root, "invalid.tsv"), writer = "unknown"),
                    error = identity)
stopifnot(inherits(invalid, "error"), !file.exists(file.path(root, "invalid.tsv")))
unlink(root, recursive = TRUE)
cat("PASS network IO utilities\n")
