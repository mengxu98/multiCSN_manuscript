suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
source("code/functions/ranking_metrics.R")

evaluate_layer_ablation <- function(raw_root = "runs/gse274113/method_fits/multiCSN", truth_root = "data/gse274113/heldout_pseudobulk_full", out = "results/gse274113/ablation") {
raw_root <- normalizePath(raw_root, mustWork = TRUE)
truth_root <- normalizePath(truth_root, mustWork = TRUE)
if (dir.exists(out)) stop("Output already exists")
dir.create(out, recursive = TRUE)
quantile_score <- function(x) rank(x, ties.method = "average") / (length(x) + 1)
project <- function(chains, column) {
  x <- copy(chains)
  x[, score := get(column)]
  x[, maximum := max(score), by = .(regulator, target)]
  x <- x[abs(score - maximum) <= 1e-12 * (1 + abs(maximum))]
  x <- x[, .(score = max(score), direction = if (uniqueN(chain_direction) == 1L) chain_direction[[1L]] else NA_real_), by = .(regulator, target)]
  x[is.finite(direction)]
}
tied_top <- function(scores, y, correct, k) {
  cutoff <- sort(scores, decreasing = TRUE)[[k]]
  above <- scores > cutoff
  equal <- scores == cutoff
  w <- as.numeric(above) + as.numeric(equal) * (k - sum(above)) / sum(equal)
  tp <- sum(w * y)
  list(precision = tp / k, recall = if (sum(y)) tp / sum(y) else NA_real_, direction_accuracy = if (tp) sum(w * correct) / tp else NA_real_, cutoff_ties = sum(equal))
}
for (day in c(7, 9, 11, 14)) {
  slug <- paste0("day_", day)
  paths <- file.path(raw_root, slug, c("tf_gene.tsv", "tf_region.tsv", "region_gene.tsv", "triplets.tsv"))
  direct <- fread(paths[[1L]])[, .(regulator, target, score = deletion_delta_bic, direction = sign(standardized_beta))]
  tr <- fread(paths[[2L]])[, .(regulator, region, tr = deletion_delta_bic)]
  rg <- fread(paths[[3L]])[, .(region, target, rg = deletion_delta_bic)]
  chains <- fread(paths[[4L]])
  tr[, tr_quantile := quantile_score(tr)]
  rg[, rg_quantile := quantile_score(rg)]
  chains <- merge(chains, tr, by = c("regulator", "region"), all.x = TRUE)
  chains <- merge(chains, rg, by = c("region", "target"), all.x = TRUE)
  chains[, `:=`(raw_min = pmin(tr, rg), quantile_min = pmin(tr_quantile, rg_quantile), geometric = sqrt(pmax(tr, 0)) * sqrt(pmax(rg, 0)))]
  variants <- list(Direct = direct, `Chain raw-min` = project(chains, "raw_min"), `Chain percentile-min` = project(chains, "quantile_min"), `Chain geometric-mean` = project(chains, "geometric"))
  for (name in names(variants)) distributions[[length(distributions) + 1L]] <- data.table(day = day, representation = name, pairs = nrow(variants[[name]]), distinct_scores = uniqueN(variants[[name]]$score), repeated_score_fraction = mean(duplicated(variants[[name]]$score)))
  d <- copy(direct); c <- copy(variants[["Chain raw-min"]])
  d[, direct_q := quantile_score(score), by = regulator]
  c[, chain_q := quantile_score(score), by = regulator]
  both <- merge(d, c, by = c("regulator", "target"), suffixes = c("_direct", "_chain"))
  both <- both[direction_direct == direction_chain, .(regulator, target, score = pmin(direct_q, chain_q), direction = direction_direct)]
  variants[["Sign-consistent intersection"]] <- both
  files <- sort(list.files(file.path(truth_root, "rna_effects"), pattern = paste0("^", slug, "__.*[.]tsv[.]gz$"), full.names = TRUE))
  for (file in files) {
    tf <- sub(paste0("^", slug, "__"), "", sub("[.]tsv[.]gz$", "", basename(file)))
    truth <- fread(file)[feature != tf & is.finite(FDR) & is.finite(logFC), .(target = feature, response = FDR <= 0.05, logFC)]
    evaluated <- lapply(variants, function(pred) {
      z <- merge(truth, pred[regulator == tf, .(target, score, direction)], by = "target", all.x = TRUE)
      z[is.na(score), `:=`(score = 0, direction = 0)]
      z[, correct := response & score > 0 & sign(logFC) == -direction]
      z
    })
    supports <- vapply(evaluated, function(x) sum(x$score > 0), integer(1))
    for (name in names(evaluated)) {
      z <- evaluated[[name]]
      ap <- threshold_average_precision(z$response, z$score)
      metrics[[length(metrics) + 1L]] <- data.table(day = day, regulator = tf, representation = name, evaluable_genes = nrow(z), responses = sum(z$response), prevalence = mean(z$response), selected = sum(z$score > 0), selected_responses = sum(z$response & z$score > 0), response_coverage = sum(z$response & z$score > 0) / sum(z$response), average_precision = ap, ap_lift = ap / mean(z$response))
      for (k in c(25L, 50L, 100L)) {
        common <- all(supports >= k)
        m <- if (common) tied_top(z$score, z$response, z$correct, k) else list(precision = NA_real_, recall = NA_real_, direction_accuracy = NA_real_, cutoff_ties = NA_integer_)
        budgets[[length(budgets) + 1L]] <- cbind(data.table(day = day, regulator = tf, representation = name, k = k, common_support_eligible = common, minimum_variant_support = min(supports)), as.data.table(m))
      }
    }
  }
}
m <- rbindlist(metrics); b <- rbindlist(budgets)
fwrite(m, file.path(out, "stratum_metrics.tsv"), sep = "\t")
fwrite(b, file.path(out, "matched_budget_metrics.tsv"), sep = "\t")
fwrite(m[, .(timepoints = .N, mean_ap = mean(average_precision), mean_ap_lift = mean(ap_lift), mean_response_coverage = mean(response_coverage), mean_selected = mean(selected)), by = .(regulator, representation)], file.path(out, "TF_macro_metrics.tsv"), sep = "\t")
fwrite(b[common_support_eligible == TRUE, .(strata = .N, TFs = uniqueN(regulator), mean_precision = mean(precision), mean_recall = mean(recall), mean_direction_accuracy = mean(direction_accuracy, na.rm = TRUE)), by = .(representation, k)], file.path(out, "matched_budget_summary.tsv"), sep = "\t")
fwrite(rbindlist(distributions), file.path(out, "score_resolution.tsv"), sep = "\t")


}

if (sys.nframe() == 0L) evaluate_layer_ablation()
