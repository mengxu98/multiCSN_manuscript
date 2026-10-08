#!/usr/bin/env Rscript

source("code/functions/ranking_metrics.R")

source("code/functions/table_io.R")

suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
result_root <- file.path(repo_root, "results", "gse274113", "method_comparison")
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)

methods <- c("multiCSN", "SCENIC+", "Pando", "LINGER")
timepoints <- c("day 7", "day 9", "day 11", "day 14")

display_method <- function(x) {
  fifelse(x == "SCENICplus", "SCENIC+", as.character(x))
}

tf_gene_formal_path <- file.path(
  repo_root, "results", "gse274113", "method_comparison", "stratum_metrics.tsv"
)
tf_gene_linger_path <- file.path(result_root, "linger_stratum_metrics.tsv")
binding_formal_path <- file.path(
  repo_root, "results", "gse274113", "multiomic_endpoint_evidence",
  "TF_region_external_binding", "TF_region_external_binding_metrics.tsv.gz"
)
binding_linger_path <- file.path(
  repo_root, "results", "gse274113", "multiomic_endpoint_evidence", "TF_region_external_binding", "LINGER", "TF_region_external_binding_metrics.tsv.gz"
)
contact_formal_path <- file.path(
  repo_root, "results", "gse274113", "multiomic_endpoint_evidence",
  "region_gene_HSPC_contact", "HSPC_contact_summary.tsv"
)
contact_linger_path <- file.path(
  repo_root, "results", "gse274113", "multiomic_endpoint_evidence", "region_gene_HSPC_contact", "LINGER", "HSPC_contact_summary.tsv"
)
chain_formal_path <- file.path(
  repo_root, "results", "gse274113", "multiomic_endpoint_evidence",
  "factorized_chain_response", "chain_response_TF_cluster_metrics.tsv"
)
chain_linger_path <- file.path(
  repo_root, "results", "gse274113", "multiomic_endpoint_evidence", "factorized_chain_response", "LINGER", "chain_response_TF_cluster_metrics.tsv"
)
input_paths <- c(
  tf_gene_formal_path, tf_gene_linger_path,
  binding_formal_path, binding_linger_path,
  contact_formal_path, contact_linger_path,
  chain_formal_path, chain_linger_path
)
if (any(!file.exists(input_paths))) {
  stop("Missing LINGER or registered comparison inputs: ", paste(input_paths[!file.exists(input_paths)], collapse = ", "))
}

tf_gene <- rbindlist(list(
  fread(tf_gene_formal_path)[
    endpoint == "TF-gene" & comparison_universe == "full_truth_feature_universe"
  ],
  fread(tf_gene_linger_path)
), use.names = TRUE, fill = TRUE)
tf_gene[, method := display_method(method)]
if (nrow(tf_gene) != 4L * 19L * 4L || !setequal(unique(tf_gene$method), methods)) {
  stop("TF-gene four-method accounting failed.")
}
tf_gene_means <- tf_gene[, .(
  response_average_precision = finite_mean(response_average_precision),
  prediction_density = finite_mean(prediction_density),
  truth_response_coverage = finite_mean(truth_response_coverage),
  finite_timepoints = sum(is.finite(response_average_precision))
), by = .(method, regulator)]
tf_gene_long <- melt(
  tf_gene_means,
  id.vars = c("method", "regulator", "finite_timepoints"),
  measure.vars = c("response_average_precision", "prediction_density"),
  variable.name = "metric", value.name = "value"
)
tf_gene_long[, `:=`(
  method = factor(method, levels = methods),
  metric = factor(
    metric,
    levels = c("response_average_precision", "prediction_density"),
    labels = c("Response AP", "Prediction density")
  )
)]
write_tsv(tf_gene_means, file.path(result_root, "TF_gene_TF_means_four_methods.tsv"))

binding <- rbindlist(list(fread(binding_formal_path), fread(binding_linger_path)), fill = TRUE)
binding[, method := display_method(method)]
common_keys <- binding[, .(
  methods_observed = uniqueN(method),
  methods_nonempty = sum(predicted_regions_in_universe > 0L)
), by = .(timepoint, regulator, evidence_source)][
  methods_observed == length(methods) & methods_nonempty == length(methods),
  .(timepoint, regulator, evidence_source)
]
binding_common <- merge(
  binding, common_keys,
  by = c("timepoint", "regulator", "evidence_source"), all = FALSE
)
binding_summary <- binding_common[, {
  total_predicted <- sum(predicted_regions_in_universe)
  total_evidence <- sum(evidence_positive_regions)
  total_true_positive <- sum(true_positive)
  list(
    common_nonempty_strata = .N,
    predicted_regions_in_universe = total_predicted,
    evidence_positive_regions = total_evidence,
    true_positive = total_true_positive,
    micro_precision = total_true_positive / total_predicted,
    micro_recall = total_true_positive / total_evidence,
    macro_enrichment_over_prevalence = finite_mean(enrichment_over_prevalence)
  )
}, by = .(method, evidence_source)]
if (!setequal(unique(binding_summary$method), methods)) {
  stop("Four-method TF-region common-nonempty accounting failed.")
}
write_tsv(
  binding_summary,
  file.path(result_root, "TF_region_binding_four_method_common_nonempty.tsv")
)

contact <- rbindlist(list(fread(contact_formal_path), fread(contact_linger_path)), fill = TRUE)
contact[, `:=`(
  target_gene_block_holm_p = p.adjust(target_gene_block_wilcoxon_p, method = "holm"),
  target_gene_block_BH_q = p.adjust(target_gene_block_wilcoxon_p, method = "BH")
)]
contact[, `:=`(
  method = factor(display_method(method), levels = methods),
  timepoint = factor(timepoint, levels = timepoints)
)]
if (nrow(contact) != 16L || anyDuplicated(contact[, .(method, timepoint)])) {
  stop("Four-method region-gene contact accounting failed.")
}
write_tsv(contact, file.path(result_root, "region_gene_contact_four_methods.tsv"))

chain <- rbindlist(list(fread(chain_formal_path), fread(chain_linger_path)), fill = TRUE)
chain[, method := display_method(method)]
if (nrow(chain) != 4L * 19L || !setequal(unique(chain$method), methods)) {
  stop("Four-method chain-response accounting failed.")
}
chain[, method := factor(method, levels = methods)]
write_tsv(chain, file.path(result_root, "chain_response_four_methods.tsv"))

summary_table <- rbindlist(list(
  tf_gene[, .(
    endpoint = "TF-gene",
    metric = "response_average_precision",
    value = finite_mean(response_average_precision),
    denominator = sum(is.finite(response_average_precision))
  ), by = method],
  tf_gene[, .(
    endpoint = "TF-gene",
    metric = "prediction_density",
    value = finite_mean(prediction_density),
    denominator = sum(is.finite(prediction_density))
  ), by = method],
  binding_summary[evidence_source == "UniBind_direct_TFBS", .(
    method, endpoint = "TF-region", metric = "UniBind_micro_precision",
    value = micro_precision, denominator = common_nonempty_strata
  )],
  binding_summary[evidence_source == "UniBind_direct_TFBS", .(
    method, endpoint = "TF-region", metric = "UniBind_micro_recall",
    value = micro_recall, denominator = common_nonempty_strata
  )],
  contact[, .(
    endpoint = "region-gene",
    metric = "mean_HSPC_contact_enrichment",
    value = mean(enrichment_over_distance_matched_expectation),
    denominator = .N
  ), by = method],
  chain[, .(
    endpoint = "TF-region-gene chain",
    metric = "mean_joint_ATAC_RNA_positive_fraction",
    value = finite_mean(joint_positive_fraction),
    denominator = sum(is.finite(joint_positive_fraction))
  ), by = method]
), use.names = TRUE)
write_tsv(summary_table, file.path(result_root, "comparison_summary.tsv"))

cat("linger_four_method_aggregation\tPASS\t", result_root, "\n", sep = "")
