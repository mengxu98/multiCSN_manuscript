#!/usr/bin/env Rscript

source("code/functions/ranking_metrics.R")

source("code/functions/table_io.R")

suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
result_root <- file.path(repo_root, "results", "gse274113", "comparison")
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)

methods <- c("multiCSN", "SCENIC+", "Pando", "LINGER")
timepoints <- c("day 7", "day 9", "day 11", "day 14")

display_method <- function(x) {
  fifelse(x == "SCENICplus", "SCENIC+", as.character(x))
}

tf_gene_formal_path <- file.path(
  repo_root, "results", "gse274113", "comparison", "stratum_metrics.tsv"
)
tf_gene_linger_path <- file.path(result_root, "linger_stratum_metrics.tsv")
binding_formal_path <- file.path(
  repo_root, "results", "gse274113", "evidence", "TF_region_external_binding_metrics.tsv.gz"
)
binding_linger_path <- file.path(
  repo_root, "results", "gse274113", "evidence", "linger_TF_region_external_binding_metrics.tsv.gz"
)
contact_formal_path <- file.path(
  repo_root, "results", "gse274113", "evidence", "HSPC_contact_summary.tsv"
)
contact_linger_path <- file.path(
  repo_root, "results", "gse274113", "evidence", "linger_HSPC_contact_summary.tsv"
)
chain_formal_path <- file.path(
  repo_root, "results", "gse274113", "evidence", "chain_response_TF_cluster_metrics.tsv"
)
chain_linger_path <- file.path(
  repo_root, "results", "gse274113", "evidence", "linger_chain_response_TF_cluster_metrics.tsv"
)

tf_gene <- rbindlist(list(
  fread(tf_gene_formal_path)[
    endpoint == "TF-gene" & comparison_universe == "full_truth_feature_universe"
  ],
  fread(tf_gene_linger_path)
), use.names = TRUE, fill = TRUE)
tf_gene[, method := display_method(method)]
truth_design <- unique(fread(file.path(result_root, "truth_stratum_summary.tsv"),
  select = c("regulator", "timepoint")))
regulators <- sort(unique(truth_design$regulator))
if (!length(regulators) || anyNA(truth_design) ||
    !fsetequal(truth_design, CJ(regulator = regulators, timepoint = timepoints))) {
  stop("Invalid input TF/timepoint truth design.")
}
if (nrow(tf_gene) != length(methods) * nrow(truth_design) ||
    anyDuplicated(tf_gene[, .(method, regulator, timepoint)]) ||
    !fsetequal(tf_gene[, .(method, regulator, timepoint)],
      CJ(method = methods, regulator = regulators, timepoint = timepoints))) {
  stop("TF-gene four-method accounting failed.")
}
tf_gene_means <- tf_gene[, .(
  response_average_precision = finite_mean(response_average_precision),
  prediction_density = finite_mean(prediction_density),
  truth_response_coverage = finite_mean(truth_response_coverage),
  finite_timepoints = sum(is.finite(response_average_precision))
), by = .(method, regulator)]
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
if (nrow(contact) != length(methods) * length(timepoints) ||
    anyDuplicated(contact[, .(method, timepoint)]) ||
    !fsetequal(contact[, .(method = as.character(method), timepoint = as.character(timepoint))],
      CJ(method = methods, timepoint = timepoints))) {
  stop("Four-method region-gene contact accounting failed.")
}
write_tsv(contact, file.path(result_root, "region_gene_contact_four_methods.tsv"))

chain <- rbindlist(list(fread(chain_formal_path), fread(chain_linger_path)), fill = TRUE)
chain[, method := display_method(method)]
if (nrow(chain) != length(methods) * length(regulators) ||
    anyDuplicated(chain[, .(method, regulator)]) ||
    !fsetequal(chain[, .(method, regulator)], CJ(method = methods, regulator = regulators))) {
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
