suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)

summarize_within_class_sensitivity <- function(run = "results/gse162170/within_class/fits", reference = ".", out = "results/gse162170/within_class") {
run <- normalizePath(run, mustWork = TRUE)
reference <- normalizePath(reference, mustWork = TRUE)
dir.create(out, recursive = TRUE)
plan <- fread(file.path(run, "task_plan.tsv"))
base <- fread(file.path(reference, "results/gse162170/network/triplets.tsv"))[, .(regulator, region, target, direction = chain_direction)]
direct <- fread(file.path(reference, "results/gse162170/network/tf_gene.tsv"))[, .(regulator, target, direction = sign(standardized_beta))]
rows <- list(); loci <- list(); retention <- list()
for (i in seq_len(nrow(plan))) {
  task <- plan[i]
  p <- file.path(run, task$task_slug)
  ch <- fread(file.path(p, "triplets.tsv.gz"))[, .(regulator, region, target, direction = chain_direction)]
  dg <- fread(file.path(p, "tf_gene.tsv.gz"))[, .(regulator, target, direction = sign(standardized_beta))]
  for (ep in c("TF-gene", "factorized-chain")) {
    x <- if (ep == "TF-gene") dg else ch
    y <- if (ep == "TF-gene") direct else base
    keys <- if (ep == "TF-gene") c("regulator", "target") else c("regulator", "region", "target")
    overlap <- merge(x, y, by = keys, suffixes = c("_within", "_pooled"))
    shared <- nrow(overlap); same <- sum(overlap$direction_within == overlap$direction_pooled)
    rows[[length(rows) + 1L]] <- data.table(cell_class = task$state, cells = task$cells, endpoint = ep, pooled_edges = nrow(y), within_class_edges = nrow(x), shared_unsigned = shared, shared_signed = same, signed_jaccard = same / (nrow(x) + nrow(y) - same), pooled_signed_recall = same / nrow(y), sign_concordance_shared = if (shared) same / shared else NA_real_)
  }
  for (gene in c("ZNF804A", "KCNC2")) {
    a <- base[target == gene]
    b <- ch[target == gene]
    ov <- merge(a, b, by = c("regulator", "region", "target"), suffixes = c("_pooled", "_within"))
    loci[[length(loci) + 1L]] <- data.table(cell_class = task$state, target = gene, pooled_chains = nrow(a), within_class_chains = nrow(b), shared_unsigned = nrow(ov), shared_signed = sum(ov$direction_pooled == ov$direction_within), pooled_signed_recall = sum(ov$direction_pooled == ov$direction_within) / nrow(a))
    b[, within_selected := TRUE]
    full <- merge(a, b, by = c("regulator", "region", "target"), all.x = TRUE, suffixes = c("_pooled", "_within"))
    full[is.na(within_selected), within_selected := FALSE]
    full[, `:=`(cell_class = task$state, cells = task$cells, retained_same_sign = within_selected & direction_pooled == direction_within)]
    retention[[length(retention) + 1L]] <- full
  }
}
fwrite(rbindlist(rows), file.path(out, "pooled_within_class_overlap.tsv"), sep = "\t")
fwrite(rbindlist(loci), file.path(out, "focal_locus_overlap.tsv"), sep = "\t")
fwrite(rbindlist(retention), file.path(out, "focal_chain_retention.tsv"), sep = "\t")

}

if (sys.nframe() == 0L) summarize_within_class_sensitivity()
