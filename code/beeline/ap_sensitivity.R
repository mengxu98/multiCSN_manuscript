suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
source("code/functions/functions.R")
root <- "results/beeline"
out <- file.path(root, "ap_sensitivity")
manifest <- fread(file.path(root, "evaluation/case_metrics.tsv"))
truths <- unlist(lapply(c("Curated", "Synthetic"), function(x) list.files(file.path("data/beeline", x), "^refNetwork.csv$", recursive = TRUE, full.names = TRUE)))
truth_map <- setNames(truths, basename(dirname(truths)))
rows <- lapply(seq_len(nrow(manifest)), function(i) {
  item <- manifest[i]
  path <- file.path(root, "networks", paste0(item$Method, "__", item$Case, ".csv"))
  truth_path <- truth_map[[item$Case]]
  metric <- fast_standard_metric_summary(as.data.frame(fread(path)), as.data.frame(fread(truth_path)))
  c(as.list(item[, .(Case, DatasetType, DatasetGroup, Method)]), metric)
})
x <- rbindlist(rows)
x[, AP_rank := frank(-AP, ties.method = "average"), by = Case]
canonical_ranks <- manifest[, .(Method, AUPRC_rank = frank(-AUPRC, ties.method = "average")), by = Case]
x[canonical_ranks, on = .(Case, Method), AUPRC_rank := i.AUPRC_rank]
summary <- rbindlist(lapply(c("all", "curated", "synthetic"), function(group) {
  z <- if (group == "all") x else x[DatasetType == group]
  z[, .(Scope = group, Cases = .N, AP = mean(AP), AUPRC = mean(AUPRC), Mean_AP_rank = mean(AP_rank), Mean_AUPRC_rank = mean(AUPRC_rank)), by = Method]
}))
groups <- x[, .(Cases = .N, AP = mean(AP), AUPRC = mean(AUPRC)), by = .(DatasetType, DatasetGroup, Method)]
paired <- merge(x[Method == "multiCSN", .(Case, multiCSN_AP = AP)], x[Method != "multiCSN", .(Case, DatasetType, Method, Comparator_AP = AP)], by = "Case")
paired[, Difference := multiCSN_AP - Comparator_AP]
paired_summary <- paired[, .(Cases = .N, Mean_difference = mean(Difference), Median_difference = median(Difference), Wins = sum(Difference > 1e-12), Ties = sum(abs(Difference) <= 1e-12), Losses = sum(Difference < -1e-12)), by = Method]
dir.create(out)
fwrite(x, file.path(out, "case_metrics.tsv"), sep = "\t")
fwrite(summary, file.path(out, "method_summary.tsv"), sep = "\t")
fwrite(groups, file.path(out, "group_summary.tsv"), sep = "\t")
fwrite(paired, file.path(out, "paired_differences.tsv"), sep = "\t")
fwrite(paired_summary, file.path(out, "paired_summary.tsv"), sep = "\t")
print(summary[Scope == "all"][order(-AP)])
