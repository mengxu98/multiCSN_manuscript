#!/usr/bin/env Rscript
prepare_complete_method_tables <- function(result_root = "results/beeline") {

suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
result_root <- normalizePath(result_root, mustWork = TRUE)
evaluation_path <- file.path(result_root, "evaluation", "case_metrics.tsv")
evaluated <- fread(evaluation_path)
if (nrow(evaluated) != 5500L || uniqueN(evaluated$Case) != 500L ||
    anyDuplicated(evaluated, by=c("Case", "Method"))) stop("Incomplete common network evaluation")
set(evaluated, j = "Dataset", value = as.character(evaluated$Case))
set(evaluated, j = "ReferenceNetwork", value = rep(NA_character_, nrow(evaluated)))
set(evaluated, j = "Species", value = rep(NA_character_, nrow(evaluated)))
output <- file.path(result_root, "metrics")
if (file.exists(output)) stop("Output already exists: ", output, call. = FALSE)
stage <- paste0(output, ".staging")
if (file.exists(stage)) stop("Staging path already exists: ", stage, call. = FALSE)
dir.create(stage, recursive = TRUE)

primary_columns <- c(
  "Dataset", "ReferenceNetwork", "Species", "Case", "DatasetGroup",
  "Method", "AUPRC", "EarlyPrecision", "EPR"
)
other_columns <- c(
  "Dataset", "ReferenceNetwork", "Species", "Case", "DatasetGroup",
  "Method", "AUROC", "ACC", "Precision", "Recall", "F1"
)
methods <- c(
  "multiCSN", "GENIE3", "GNIPLR", "GRNBoost2", "LEAP", "PIDC",
  "PPCOR", "SCODE", "SCORPION", "SCSGL", "SINCERITIES"
)
primary_groups <- list()
other_groups <- list()

for (dataset_type in c("curated", "synthetic")) {
  selected <- evaluated[DatasetType == dataset_type]
  primary <- selected[, ..primary_columns]
  other <- selected[, ..other_columns]
  setorder(primary, Case, Method)
  setorder(other, Case, Method)
  expected_cases <- if (dataset_type == "curated") 200L else 300L
  if (nrow(primary) != expected_cases * 11L || nrow(other) != nrow(primary) ||
    uniqueN(primary$Case) != expected_cases ||
    !setequal(unique(primary$Method), methods) ||
    anyDuplicated(primary, by = c("Case", "Method")) ||
    anyDuplicated(other, by = c("Case", "Method"))) {
    stop("Integrated BEELINE metric contract failed for ", dataset_type, ".")
  }
  fwrite(primary, file.path(stage, paste0(dataset_type, "__primary_metrics.csv")))
  fwrite(other, file.path(stage, paste0(dataset_type, "__other_metrics.csv")))
  primary_groups[[dataset_type]] <- copy(primary)[, DatasetType := dataset_type]
  other_groups[[dataset_type]] <- copy(other)[, DatasetType := dataset_type]
}

primary_combined <- rbindlist(primary_groups, use.names = TRUE)
other_combined <- rbindlist(other_groups, use.names = TRUE)
setorder(primary_combined, DatasetType, Case, Method)
setorder(other_combined, DatasetType, Case, Method)
fwrite(other_combined, file.path(stage, "combined__other_metrics.csv"))

panel_a <- copy(primary_combined)
fwrite(
  panel_a, file.path(stage, "case_metrics.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
ranks <- melt(
  panel_a,
  id.vars = c("Dataset", "Case", "DatasetGroup", "DatasetType", "Method"),
  measure.vars = c("AUPRC", "EarlyPrecision", "EPR"),
  variable.name = "Metric", value.name = "Value"
)
if (any(!is.finite(ranks$Value))) stop("Nonfinite primary metrics require review before ranking")
ranks[, Rank := frank(-Value, ties.method = "average"), by = .(Case, Metric)]
setorder(ranks, Metric, Case, Rank, Method)
rank_summary <- ranks[, .(
  n_cases = .N, mean_rank = mean(Rank), sd_rank = sd(Rank),
  median_rank = median(Rank), q25_rank = quantile(Rank, 0.25),
  q75_rank = quantile(Rank, 0.75)
), by = .(Method, Metric)]
setorder(rank_summary, Metric, mean_rank, Method)
fwrite(ranks, file.path(stage, "case_ranks.tsv"), sep = "\t", quote = FALSE)
fwrite(rank_summary, file.path(stage, "rank_summary.tsv"), sep = "\t", quote = FALSE)

# Source Data uses these two summaries; derive them in the same pass as ranks.
method_summary <- merge(
  primary_combined[, .(
    cases = .N, mean_auprc = mean(AUPRC),
    curated_mean_auprc = mean(AUPRC[DatasetType == "curated"]),
    synthetic_mean_auprc = mean(AUPRC[DatasetType == "synthetic"])
  ), by = .(method = Method)],
  rank_summary[Metric == "AUPRC", .(method = Method, mean_within_case_auprc_rank = mean_rank)],
  by = "method"
)
setorder(method_summary, mean_within_case_auprc_rank, method)
fwrite(method_summary, file.path(result_root, "evaluation", "beeline_method_summary.tsv"),
       sep = "\t", quote = FALSE)
case_rows <- lapply(sort(unique(evaluated$Case)), function(case) {
  row <- evaluated[Case == case][1L]
  directory <- if (row$DatasetType == "curated") "Curated" else "Synthetic"
  truth <- fread(file.path(repo, "data/beeline", directory, row$DatasetGroup, case, "refNetwork.csv"))
  names(truth)[1:2] <- c("regulator", "target")
  truth <- unique(truth[regulator != target, .(regulator, target)])
  genes <- unique(c(truth$regulator, truth$target))
  model <- sub("-[0-9]+-[0-9]+(-[0-9]+)?$", "", case)
  fields <- as.integer(strsplit(substring(case, nchar(model) + 2L), "-", fixed = TRUE)[[1L]])
  if (!length(fields) %in% c(2L, 3L) || anyNA(fields) || nrow(truth) != row$GoldEdges) {
    stop("Case/truth identity changed: ", case)
  }
  data.table(case_id = case, model = model, dataset_type = row$DatasetType,
             cell_count = fields[[1L]], replicate = fields[[2L]],
             dropout_percent = if (length(fields) == 3L) fields[[3L]] else 0L,
             gene_count = length(genes), truth_edge_count = nrow(truth),
             directed_nonself_universe = length(genes) * (length(genes) - 1L),
             gene_count_basis = "unique gene identifiers in directed nonself truth edges")
})
fwrite(rbindlist(case_rows), file.path(result_root, "evaluation", "beeline_case_manifest.tsv"),
       sep = "\t", quote = FALSE)

if (!file.rename(stage, output)) stop("Unable to promote BEELINE metrics.")
cat("beeline_metrics\tPASS\t", output, "\n", sep = "")
}



if (sys.nframe() == 0L) prepare_complete_method_tables()
