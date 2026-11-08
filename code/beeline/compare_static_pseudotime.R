args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply a new output directory")
output <- args[[1L]]
if (file.exists(output)) stop("Output already exists: ", output)
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
source("code/functions/functions.R")
source("code/methods/multiCSN.R")
dir.create(output, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
writeLines(find.package("inferCSN"), file.path(output, "package_path.txt"))
baseline <- fread("results/beeline/ap_sensitivity/case_metrics.tsv")
keys <- c("Case", "DatasetType", "DatasetGroup")
paths <- sort(unlist(lapply(c("Curated", "Synthetic"), function(group) {
  list.files(file.path("data/beeline", group), "^ExpressionData[.]csv$",
             recursive = TRUE, full.names = TRUE)
})))
cases <- rbindlist(lapply(paths, function(path) {
  data.table(Case = basename(dirname(path)),
             DatasetType = tolower(basename(dirname(dirname(dirname(path))))),
             DatasetGroup = basename(dirname(dirname(path))), expression = path)
}))
stopifnot(nrow(cases) == 500L, !anyDuplicated(cases$Case),
          sum(cases$DatasetType == "curated") == 200L,
          sum(cases$DatasetType == "synthetic") == 300L,
          nrow(baseline) == 5500L, !anyDuplicated(baseline, by = c("Case", "Method")),
          fsetequal(cases[, ..keys], unique(baseline[, ..keys])))
fwrite(cases, file.path(output, "cases.tsv"), sep = "\t")
metrics <- c("AP", "AUPRC", "AUROC", "EarlyPrecision", "EPR", "ACC", "Precision", "Recall", "F1")
runtime <- list()
results <- list()
checks <- list()
canonical <- function(x) {
  x <- as.data.table(x)
  setorder(x, regulator, target)
  x
}
for (repeat_id in seq_len(5L)) {
  for (case_id in seq_len(nrow(cases))) {
    item <- cases[case_id]
    directory <- dirname(item$expression)
    load_start <- Sys.time()
    expression <- t(as.matrix(read.csv(item$expression, row.names = 1L, check.names = FALSE)))
    pseudotime <- read.csv(file.path(directory, "PseudoTime.csv"), row.names = 1L, check.names = FALSE)
    truth <- as.data.frame(fread(file.path(directory, "refNetwork.csv")))
    load_seconds <- as.numeric(difftime(Sys.time(), load_start, units = "secs"))
    stopifnot(nrow(expression) == nrow(pseudotime), identical(rownames(expression), rownames(pseudotime)))
    modes <- c("static", "pseudotime")
    if ((repeat_id + case_id) %% 2L == 0L) modes <- rev(modes)
    for (mode in modes) {
      destination <- file.path(output, paste0("repeat_", repeat_id), mode,
                               paste0(item$Case, ".csv"))
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      set.seed(1L)
      started <- Sys.time()
      network <- run_multiCSN(expression, if (mode == "static") NULL else pseudotime,
                             cores = 1L, filename = destination)
      seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      stopifnot(is.data.frame(network), identical(names(network), c("regulator", "target", "weight")),
                nrow(network) > 0L, all(is.finite(network$weight)),
                all(network$regulator != network$target),
                !anyDuplicated(network[c("regulator", "target")]),
                all(network$regulator %in% colnames(expression)),
                all(network$target %in% colnames(expression)), is.finite(seconds), seconds > 0)
      saved <- fread(destination)
      stopifnot(isTRUE(all.equal(canonical(network), canonical(saved), tolerance = 1e-12)))
      if (repeat_id > 1L) {
        first <- fread(file.path(output, "repeat_1", mode, paste0(item$Case, ".csv")))
        stopifnot(isTRUE(all.equal(canonical(saved), canonical(first), tolerance = 1e-12)))
      }
      runtime[[length(runtime) + 1L]] <- cbind(item[, ..keys], data.table(
        Mode = mode, Repeat = repeat_id, Order = match(mode, modes),
        RuntimeSeconds = seconds, LoadSeconds = load_seconds,
        SelectedEdges = nrow(saved), TimingScope = "inference_and_export", Status = "complete"))
      if (repeat_id == 1L) {
        score <- fast_standard_metric_summary(as.data.frame(saved), truth)
        stopifnot(all(is.finite(unlist(score[metrics]))))
        results[[length(results) + 1L]] <- c(as.list(item[, ..keys]),
          list(Method = paste0("multiCSN (", mode, ")")), score)
        if (mode == "pseudotime") {
          historical <- fread(file.path("results/beeline/networks", paste0("multiCSN__", item$Case, ".csv")))
          historical_score <- fast_standard_metric_summary(as.data.frame(historical), truth)
          recorded <- baseline[Case == item$Case & Method == "multiCSN"]
          score_error <- max(abs(unlist(score[metrics]) - unlist(historical_score[metrics])))
          recorded_error <- max(abs(unlist(historical_score[metrics]) - as.numeric(recorded[, ..metrics])))
          identical_network <- isTRUE(all.equal(canonical(saved), canonical(historical), tolerance = 1e-10))
          checks[[length(checks) + 1L]] <- data.table(Case = item$Case,
            NetworkMatch = identical_network, ReplayMetricMaxError = score_error,
            BaselineMetricMaxError = recorded_error)
          stopifnot(identical_network, score_error < 1e-10, recorded_error < 1e-10)
        }
      }
    }
    if (case_id %% 50L == 0L) {
      fwrite(rbindlist(runtime), file.path(output, "case_runtime.tsv"), sep = "\t")
      message("repeat ", repeat_id, "/5 cases ", case_id, "/500")
    }
  }
}
x <- rbindlist(results, fill = TRUE)
runtimes <- rbindlist(runtime)
stopifnot(nrow(x) == 1000L, nrow(runtimes) == 5000L,
          all(runtimes[, .N, by = .(Mode, Repeat)]$N == 500L))
fwrite(x, file.path(output, "mode_case_metrics.tsv"), sep = "\t")
fwrite(rbindlist(checks), file.path(output, "baseline_replay_checks.tsv"), sep = "\t")
fwrite(runtimes, file.path(output, "case_runtime.tsv"), sep = "\t")
combined <- rbindlist(list(baseline[Method != "multiCSN"], x), use.names = TRUE, fill = TRUE)
for (metric in metrics) {
  combined[, (paste0(metric, "_rank12")) := frank(-get(metric), ties.method = "average"), by = Case]
}
fwrite(combined, file.path(output, "all_method_case_metrics.tsv"), sep = "\t")
summary <- rbindlist(lapply(c("all", "curated", "synthetic"), function(scope) {
  z <- if (scope == "all") combined else combined[DatasetType == scope]
  z[, c(list(Scope = scope, Cases = .N), lapply(.SD, mean)),
    by = Method, .SDcols = c(metrics, paste0(metrics, "_rank12"))]
}))
fwrite(summary, file.path(output, "method_summary.tsv"), sep = "\t")
group_summary <- combined[, c(list(Cases = .N), lapply(.SD, mean)),
  by = .(DatasetType, DatasetGroup, Method), .SDcols = metrics]
fwrite(group_summary, file.path(output, "group_summary.tsv"), sep = "\t")
static <- combined[Method == "multiCSN (static)"]
others <- combined[Method != "multiCSN (static)"]
paired <- merge(melt(static, id.vars = keys, measure.vars = metrics,
                     variable.name = "Metric", value.name = "Static"),
                melt(others, id.vars = c(keys, "Method"), measure.vars = metrics,
                     variable.name = "Metric", value.name = "Comparator"), by = c(keys, "Metric"))
paired[, Difference := Static - Comparator]
fwrite(paired, file.path(output, "paired_differences.tsv"), sep = "\t")
paired_summary <- rbindlist(lapply(c("all", "curated", "synthetic"), function(scope) {
  z <- if (scope == "all") paired else paired[DatasetType == scope]
  z[, .(Scope = scope, Cases = .N, MeanDifference = mean(Difference),
        MedianDifference = median(Difference), Wins = sum(Difference > 1e-12),
        Ties = sum(abs(Difference) <= 1e-12), Losses = sum(Difference < -1e-12)),
    by = .(Method, Metric)]
}))
fwrite(paired_summary, file.path(output, "paired_summary.tsv"), sep = "\t")
totals <- runtimes[, .(Cases = .N, TotalSeconds = sum(RuntimeSeconds),
                       MedianCaseSeconds = median(RuntimeSeconds)), by = .(Mode, Repeat)]
fwrite(totals, file.path(output, "runtime_repeats.tsv"), sep = "\t")
fwrite(totals[, .(Repeats = .N, MedianTotalSeconds = median(TotalSeconds),
                 MinimumTotalSeconds = min(TotalSeconds), MaximumTotalSeconds = max(TotalSeconds)), by = Mode],
       file.path(output, "runtime_summary.tsv"), sep = "\t")
profile_dir <- file.path(output, "profile")
dir.create(profile_dir)
for (mode in c("static", "pseudotime")) {
  Rprof(file.path(profile_dir, paste0(mode, ".Rprof")), interval = 0.001)
  for (path in paths[seq_len(50L)]) {
    expression <- t(as.matrix(read.csv(path, row.names = 1L, check.names = FALSE)))
    pseudotime <- read.csv(file.path(dirname(path), "PseudoTime.csv"), row.names = 1L, check.names = FALSE)
    invisible(run_multiCSN(expression, if (mode == "static") NULL else pseudotime,
                           cores = 1L, filename = file.path(profile_dir, paste0(mode, ".csv"))))
  }
  Rprof(NULL)
  writeLines(capture.output(summaryRprof(file.path(profile_dir, paste0(mode, ".Rprof")))),
             file.path(profile_dir, paste0(mode, "_summary.txt")))
}
writeLines("VERIFIED: 500 static cases, 500 matching historical pseudotime replays, five paired timing repeats per mode.",
           file.path(output, "COMPLETE.txt"))
print(summary[Scope == "all"][order(-AP), .(Method, AP, AUPRC, AP_rank12, AUPRC_rank12)])
