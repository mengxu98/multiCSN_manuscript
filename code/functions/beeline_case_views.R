# Shared by manuscript plots and the Python Source Data exporter.
# Reads formal saved networks; never reads or updates topology snapshots.
read_beeline_case <- function(repo, case) {
  methods <- c("GENIE3", "GNIPLR", "GRNBoost2", "LEAP", "PIDC", "PPCOR",
               "SCODE", "SCORPION", "SCSGL", "SINCERITIES", "multiCSN")
  group <- sub("-.*", "", case)
  truth <- data.table::fread(file.path(repo, "data/beeline/Curated", group, case, "refNetwork.csv"))[, 1:3]
  data.table::setnames(truth, c("regulator", "target", "sign"))
  truth <- unique(truth[regulator != target], by = c("regulator", "target"))
  truth[, sign := ifelse(sign == "+", 1L, -1L)]
  metrics <- data.table::fread(file.path(repo, "results/beeline/evaluation/case_metrics.tsv"))
  metrics <- metrics[Case == case & DatasetType == "curated" & DatasetGroup == group]
  edges <- summaries <- vector("list", length(methods))
  for (i in seq_along(methods)) {
    method <- methods[[i]]
    n <- data.table::fread(file.path(repo, "results/beeline/networks", paste0(method, "__", case, ".csv")))[, 1:3]
    data.table::setnames(n, c("regulator", "target", "raw_weight"))
    n <- n[is.finite(raw_weight) & raw_weight != 0 & regulator != target]
    n[, magnitude := abs(raw_weight)]
    data.table::setorder(n, -magnitude, regulator, target)
    n <- unique(n, by = c("regulator", "target"))
    n[, within_network_signed_rank := sign(raw_weight) * rank(magnitude, ties.method = "average") / .N]
    if (case == "HSC-2000-3-30") n <- head(n, 33L)
    n[, Method := method]
    edges[[i]] <- n[, .(Method, regulator, target, raw_weight, within_network_signed_rank)]
    hits <- merge(n, truth, by = c("regulator", "target"))
    tp <- nrow(hits)
    precision <- if (nrow(n)) tp / nrow(n) else 0
    recall <- tp / nrow(truth)
    metric <- metrics[Method == method]
    summaries[[i]] <- data.table::data.table(
      Method = method, predicted_edges = nrow(n), truth_edges_recovered = tp,
      positive_edges = sum(n$raw_weight > 0), negative_edges = sum(n$raw_weight < 0),
      display_transform = "signed ordinal rank within this network only",
      AUPRC = metric$AUPRC, EarlyPrecision = metric$EarlyPrecision, EPR = metric$EPR,
      precision = precision, recall = recall,
      F1 = if (precision + recall > 0) 2 * precision * recall / (precision + recall) else 0,
      sign_concordant = sum(sign(hits$raw_weight) == hits$sign))
  }
  list(truth = truth, edges = data.table::rbindlist(edges), summary = data.table::rbindlist(summaries))
}

write_beeline_case <- function(repo, case, output) {
  suppressPackageStartupMessages(library(data.table))
  view <- read_beeline_case(repo, case)
  for (name in names(view)) data.table::fwrite(view[[name]], file.path(output, paste0(name, ".tsv")), sep = "\t")
}
