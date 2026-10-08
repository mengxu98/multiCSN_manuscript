source("code/functions/ranking_metrics.R")

# Independent expectations for the shared perturbation ranking metrics.
truth <- c(TRUE, FALSE, TRUE)
scores <- c(3, 2, 1)
stopifnot(abs(threshold_average_precision(truth, scores) - 5 / 6) < 1e-12,
          tie_safe_auroc(truth, scores) == 0.5,
          expected_precision_at_k(truth, scores, 2L) == 0.5,
          abs(expected_precision_at_k(truth, rep(1, 3), 1L) - 2 / 3) < 1e-12,
          tie_safe_auroc(truth, rep(1, 3)) == 0.5,
          is.na(threshold_average_precision(rep(FALSE, 3), scores)),
          is.na(tie_safe_auroc(rep(TRUE, 3), scores)),
          finite_mean(c(1, 3, NA_real_, Inf)) == 2,
          is.na(finite_mean(numeric())))

# Products of class counts can exceed R's 32-bit integer range.
large_truth <- rep(c(TRUE, FALSE), each = 50000L)
large_scores <- rep(c(1, 0), each = 50000L)
stopifnot(tie_safe_auroc(large_truth, large_scores) == 1,
          tie_safe_auroc(large_truth, rep(0, length(large_truth))) == 0.5,
          ranking_metrics_tie_aware(large_scores, large_truth)$auroc == 1)

ordinal <- signed_network_ordinal_weight(c(4, -2, 2 + 5e-13, 0, NA_real_))
stopifnot(
  identical(sign(ordinal[1:3]), c(1, -1, 1)),
  identical(ordinal[2], -ordinal[3]),
  ordinal[[1]] > abs(ordinal[[2]]),
  identical(ordinal[4:5], c(0, 0))
)

scores_a <- c(1, 0, 0, 0, 0)
labels_a <- c(TRUE, TRUE, FALSE, TRUE, FALSE)
scores_b <- c(1, 0, 0, 0, 0)
labels_b <- c(TRUE, FALSE, TRUE, FALSE, TRUE)

metrics_a <- ranking_metrics_tie_aware(scores_a, labels_a, top_k = 2L)
metrics_b <- ranking_metrics_tie_aware(scores_b, labels_b, top_k = 2L)

stopifnot(isTRUE(all.equal(metrics_a$average_precision, metrics_b$average_precision)))
stopifnot(isTRUE(all.equal(metrics_a$trapezoidal_auprc, metrics_b$trapezoidal_auprc)))
stopifnot(isTRUE(all.equal(metrics_a$auroc, metrics_b$auroc)))
stopifnot(metrics_a$topk$n_selected[[1]] == 2L)
stopifnot(isTRUE(all.equal(metrics_a$topk$precision[[1]], 0.75)))

contract_metrics <- ranking_metrics_tie_aware(
  c(0.9, 0.7, 0.7, 0.2, 0, 0),
  c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE)
)
stopifnot(isTRUE(all.equal(contract_metrics$trapezoidal_auprc, 61 / 72, tolerance = 1e-14)))
stopifnot(isTRUE(all.equal(contract_metrics$auroc, 5 / 6, tolerance = 1e-14)))

empty_metrics <- ranking_metrics_tie_aware(rep(0, 6), c(TRUE, TRUE, rep(FALSE, 4)))
stopifnot(isTRUE(all.equal(empty_metrics$average_precision, 2 / 6, tolerance = 1e-14)))
stopifnot(isTRUE(all.equal(empty_metrics$trapezoidal_auprc, (1 + 2 / 6) / 2, tolerance = 1e-14)))

cat("ranking metric regression tests passed\n")
