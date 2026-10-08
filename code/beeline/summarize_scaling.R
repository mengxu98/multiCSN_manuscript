#!/usr/bin/env Rscript
summarize_scaling <- function(input = "results/beeline/core_scaling", output = input) {
x <- read.csv(file.path(input, "core_scaling_runs.csv"))
cores <- c(1L, 2L, 4L, 6L, 8L, 10L)
stopifnot(nrow(x) == 60L, length(unique(x$Case)) == 10L,
          !anyNA(x), !anyDuplicated(x[c("Case", "Cores")]),
          all(is.finite(x$RuntimeSeconds)), all(x$RuntimeSeconds > 0),
          all(is.finite(x$MemoryTreePeakMb)), all(x$MemoryTreePeakMb > 0),
          all(x$TimingScope == "public_inference_only"))
groups <- split(x, x$Case)
for (g in groups) stopifnot(identical(sort(g$Cores), cores),
                           length(unique(g$Scenario)) == 1L,
                           length(unique(g$AxisValue)) == 1L)
cols <- c("Case", "Scenario", "AxisValue", "RuntimeSeconds", "MemoryTreePeakMb", "OutputEdges")
comparison <- merge(x[x$Cores == 1L, cols], x[x$Cores == 10L, cols],
                    by = c("Case", "Scenario", "AxisValue"), suffixes = c("_1", "_10"))
comparison <- comparison[c("Case", "Scenario", "AxisValue", "RuntimeSeconds_1", "RuntimeSeconds_10",
                           "MemoryTreePeakMb_1", "MemoryTreePeakMb_10", "OutputEdges_1", "OutputEdges_10")]
comparison$Speedup_1_to_10 <- comparison$RuntimeSeconds_1 / comparison$RuntimeSeconds_10
# Full sampled peak RSS, without subtracting MemoryTreeBeforeMb.
comparison$RssRatio_10_to_1 <- comparison$MemoryTreePeakMb_10 / comparison$MemoryTreePeakMb_1
comparison <- comparison[order(comparison$Case, method = "radix"), ]
summary <- data.frame(
  statistic = c("complete_runs", "complete_cases", "median_speedup_1_to_10",
                "minimum_speedup_1_to_10", "maximum_speedup_1_to_10",
                "median_rss_ratio_10_to_1", "minimum_rss_ratio_10_to_1", "maximum_rss_ratio_10_to_1"),
  value = c(nrow(x), length(groups), median(comparison$Speedup_1_to_10),
            min(comparison$Speedup_1_to_10), max(comparison$Speedup_1_to_10),
            median(comparison$RssRatio_10_to_1), min(comparison$RssRatio_10_to_1), max(comparison$RssRatio_10_to_1)))
invariance <- do.call(rbind, lapply(groups, function(g) {
  files <- file.path(input, "runs", g$Case[[1L]], paste0("cores_", cores), "network.rds")
  exact <- NA
  if (all(file.exists(files))) {
    networks <- lapply(files, readRDS)
    stopifnot(all(vapply(networks, nrow, integer(1L)) == g$OutputEdges[match(cores, g$Cores)]))
    exact <- all(vapply(networks, identical, logical(1L), networks[[1L]]))
  }
  counts <- unique(g$OutputEdges)
  data.frame(Case = g$Case[[1L]], CoreSettings = nrow(g), EdgeCounts = length(counts),
             OutputEdges = if (length(counts) == 1L) counts else NA, FullOutputExact = exact)
}))
invariance <- invariance[order(invariance$Case, method = "radix"), ]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(comparison, file.path(output, "one_to_ten_core_comparison.csv"), row.names = FALSE)
write.csv(summary, file.path(output, "scaling_summary_statistics.csv"), row.names = FALSE)
write.csv(invariance, file.path(output, "edge_count_invariance.csv"), row.names = FALSE, na = "NA")
message("Summarized 60 runs; full-output equality verified for ", sum(!is.na(invariance$FullOutputExact)), "/10 cases.")
}

if (sys.nframe() == 0L) summarize_scaling()
