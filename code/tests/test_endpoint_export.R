# Synthetic adapter-contract fixture, not a biological rerun or historical replay.
source("code/gse274113/export_endpoints.R")
root <- tempfile("endpoint-export-")
dir.create(root)
put <- function(x, directory, name) write_tsv(as.data.table(x), file.path(directory, name))
for (method in c("multiCSN", "Pando", "SCENICplus")) {
  for (slug in c("day_7", "day_9", "day_11", "day_14")) {
    directory <- file.path(root, "input", method, slug)
    dir.create(directory, recursive = TRUE)
    if (method == "multiCSN") {
      edge <- data.table(regulator = "TF", region = "chr1:1-10", target = "G", weight = -0.75)
      for (file in c("tf_gene.tsv", "tf_region.tsv", "region_gene.tsv", "triplets.tsv"))
        put(edge, directory, file)
      put(data.table(target = c("G", "OTHER"), status = "ok"), directory, "tf_gene_accounting.tsv")
    } else if (method == "Pando") {
      edge <- data.table(tf = "TF", region = "chr1:1-10", target = "G", estimate = -2, padj = 0.01)
      put(edge, directory, "native_module_metadata.tsv")
      put(edge, directory, "retained_native_triplets.tsv.gz")
      put(edge[, .(tf, region)], directory, "retained_tf_region_binary.tsv.gz")
      put(edge[, .(region, target)], directory, "retained_region_gene_binary.tsv.gz")
      put(data.table(gene = c("G", "OTHER")), directory, "native_variable_features.tsv")
    } else {
      edge <- data.table(TF = "TF", Region = "chr1:1-10", Gene = "G",
        importance_R2G = 4, rho_R2G = 0.2, importance_TF2G = 3, rho_TF2G = -0.25, triplet_rank = 7)
      put(edge, directory, "eRegulon_direct.tsv")
      put(data.table(target = c("G", "OTHER")), directory, "tf_to_gene_adj.tsv")
    }
  }
  destination <- file.path(root, "canonical", method)
  export_endpoints(method, file.path(root, "input", method), destination)
  expected <- switch(method, multiCSN = 0.75, Pando = 2, SCENICplus = 3)
  for (slug in c("day_7", "day_9", "day_11", "day_14")) {
    tg <- read_fit_checkpoint(file.path(destination, slug, "tf_gene.tsv.gz"))
    stopifnot(tg$rank_magnitude == expected, tg$direction_sign == -1,
              nrow(read_fit_checkpoint(file.path(destination, slug, "tf_gene_target_eligibility.tsv.gz"))) == 2)
    tr <- read_fit_checkpoint(file.path(destination, slug, "tf_region.tsv.gz"))
    if (method != "multiCSN") stopifnot(!tr$score_available, tr$rank_magnitude == 0, is.na(tr$direction_sign))
  }
  failed <- tryCatch({export_endpoints(method, file.path(root, "input", method), destination); FALSE},
                     error = function(e) grepl("Output already exists", conditionMessage(e)))
  stopifnot(failed)
}
conflict <- data.table(tf = c("TF", "TF"), gene = "G", score = 2, direction = c(-1, 1))
stopifnot(inherits(tryCatch(collapse_best(conflict, c("tf", "gene"), "score", "direction", "fixture"),
                           error = identity), "error"))
unlink(root, recursive = TRUE)
cat("PASS three-method/four-timepoint endpoint exports, eligibility, unsigned relations and overwrite/conflict rejection\n")
