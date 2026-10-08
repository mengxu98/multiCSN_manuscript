suppressPackageStartupMessages(library(data.table))
source("code/functions/network_io.R")
# Exercise the real loader without starting a full perturbation evaluation.
expressions <- parse("code/gse274113/evaluate_tf_gene_comparators.R")
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("load_canonical"))) eval(expr)
}
count <- 0L
for (method in c("multiCSN", "Pando", "SCENICplus")) {
  for (slug in c("day_7", "day_9", "day_11", "day_14")) {
    root <- file.path("results/gse274113/networks", method)
    result <- load_canonical(root, slug)
    original <- read_fit_checkpoint(file.path(root, slug, "tf_gene.tsv.gz"))
    stopifnot(identical(result$table$rank_magnitude, original$rank_magnitude),
              identical(result$table$direction_sign, original$direction_sign),
              identical(result$table$score_available, original$score_available),
              nrow(result$table) == nrow(original))
    count <- count + 1L
  }
}
cat("PASS", count, "real canonical TF-gene networks loaded without score transformation\n")
