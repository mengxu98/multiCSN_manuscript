#!/usr/bin/env Rscript
# Primary paired LINGER-vs-multiCSN comparison for the held-out perturbation
# metrics quoted in the manuscript. Method (see Methods): per TF, average the
# finite paired differences across time points, then run a two-sided one-sample
# signed-rank test on the TF-level differences (asymptotic), Holm-corrected
# across Pando, SCENIC+ and LINGER separately within each metric, matching Fig. 3.
suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)

multicsn_path <- "results/gse274113/method_comparison/stratum_metrics.tsv"
linger_path <- "results/gse274113/method_comparison/linger_stratum_metrics.tsv"
output_path <- "results/gse274113/method_comparison/linger_vs_multicsn_paired_tests.tsv"

metrics <- c(response_average_precision = "Average precision",
             truth_response_coverage = "Response coverage",
             response_precision_at_truth_k = "Top-k precision")
multicsn <- fread(multicsn_path)[endpoint == "TF-gene" & comparison_universe == "full_truth_feature_universe"]
linger <- fread(linger_path)[endpoint == "TF-gene" & comparison_universe == "full_truth_feature_universe"]
multi_only <- multicsn[method == "multiCSN"]
rows <- rbindlist(lapply(names(metrics), function(metric) {
  a <- multi_only[, .(regulator, timepoint, value = as.numeric(get(metric)))]
  b <- linger[, .(regulator, timepoint, value = as.numeric(get(metric)))]
  merged <- merge(a, b, by = c("regulator", "timepoint"), suffixes = c("_multicsn", "_linger"))
  merged[, difference := value_multicsn - value_linger]
  cluster <- merged[is.finite(difference), .(mean_difference = mean(difference)), by = regulator]
  d <- cluster$mean_difference
  p <- if (length(d) >= 3L && any(d != 0)) {
    wilcox.test(d, mu = 0, paired = FALSE, exact = FALSE, conf.int = FALSE)$p.value
  } else NA_real_
  data.table(comparator = "LINGER", reference = "multiCSN", endpoint = "TF-gene",
             comparison_universe = "full_truth_feature_universe", metric = metric,
             metric_label = metrics[[metric]],
             difference_definition = "multiCSN minus LINGER",
             regulator_clusters_finite = length(d),
             mean_cluster_difference = mean(d), median_cluster_difference = median(d),
             wilcoxon_signed_rank_p = p)
}))
other <- fread("results/gse274113/method_cluster_inference/paired_method_cluster_tests.tsv")
other <- other[reference == "multiCSN" & endpoint == "TF-gene" &
  comparison_universe == "full_truth_feature_universe" &
  comparator %in% c("Pando", "SCENIC+") & metric %in% names(metrics),
  .(metric, comparator, wilcoxon_signed_rank_p)]
family <- rbindlist(list(other, rows[, .(metric, comparator, wilcoxon_signed_rank_p)]))
stopifnot(nrow(family) == 9L, all(family[, .N, by = metric]$N == 3L))
family[, holm_p_three_comparators := p.adjust(wilcoxon_signed_rank_p, method = "holm"), by = metric]
rows <- merge(rows, family[, .(metric, comparator, holm_p_three_comparators)],
              by = c("metric", "comparator"), sort = FALSE)
rows[, correction_family := "Pando, SCENIC+, LINGER within each metric"]
fwrite(rows, output_path, sep = "\t", quote = FALSE, na = "NA")
print(rows, digits = 4)
