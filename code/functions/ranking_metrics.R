signed_network_ordinal_weight <- function(values, tolerance = 1e-12) {
  values <- suppressWarnings(as.numeric(values))
  out <- numeric(length(values))
  keep <- is.finite(values) & values != 0
  if (!any(keep)) {
    return(out)
  }

  kept_index <- which(keep)
  magnitude <- abs(values[kept_index])
  ord <- order(-magnitude, kept_index)
  sorted_magnitude <- magnitude[ord]
  group <- integer(length(ord))
  group[[1]] <- 1L
  anchor <- sorted_magnitude[[1]]
  if (length(ord) > 1L) {
    for (i in 2:length(ord)) {
      if (abs(sorted_magnitude[[i]] - anchor) >= tolerance) {
        group[[i]] <- group[[i - 1L]] + 1L
        anchor <- sorted_magnitude[[i]]
      } else {
        group[[i]] <- group[[i - 1L]]
      }
    }
  }
  n_groups <- max(group)
  out[kept_index[ord]] <- sign(values[kept_index[ord]]) *
    (1 - (group - 0.5) / n_groups)
  out
}

ranking_metrics_tie_aware <- function(scores, labels, top_k = c(25L, 50L, 100L, 200L, 500L, 1000L)) {
  scores <- suppressWarnings(as.numeric(scores))
  labels <- as.logical(labels)
  if (length(scores) != length(labels)) {
    stop("scores and labels must have the same length.", call. = FALSE)
  }
  if (length(scores) == 0L || any(!is.finite(scores)) || anyNA(labels)) {
    stop("Scores must be finite and labels must be non-missing.", call. = FALSE)
  }

  n_positive <- sum(labels)
  n_negative <- length(labels) - n_positive
  if (n_positive == 0L || n_negative == 0L) {
    return(list(
      auroc = NA_real_,
      average_precision = NA_real_,
      trapezoidal_auprc = NA_real_,
      topk = data.frame(
        top_k_requested = as.integer(top_k),
        n_selected = NA_integer_,
        cutoff_score = NA_real_,
        precision = NA_real_,
        recall = NA_real_,
        f1 = NA_real_
      )
    ))
  }

  ord <- order(-scores, seq_along(scores))
  ordered_scores <- scores[ord]
  ordered_labels <- labels[ord]
  block_end <- c(which(diff(ordered_scores) != 0), length(ordered_scores))
  block_start <- c(1L, head(block_end, -1L) + 1L)
  n_selected <- block_end
  true_selected <- cumsum(ordered_labels)[block_end]
  recall <- true_selected / n_positive
  precision <- true_selected / n_selected

  average_precision <- sum(c(recall[1L], diff(recall)) * precision)

  pr_recall <- c(0, recall)
  pr_precision <- c(1, precision)
  trapezoidal_auprc <- sum(
    diff(pr_recall) *
      (head(pr_precision, -1L) + tail(pr_precision, -1L)) / 2
  )

  ranks <- rank(scores, ties.method = "average")
  auroc <- (sum(ranks[labels]) - n_positive * (n_positive + 1) / 2) /
    (as.double(n_positive) * n_negative)

  topk <- lapply(as.integer(top_k), function(k) {
    k <- min(k, length(scores))
    cutoff_score <- ordered_scores[[k]]
    above <- ordered_scores > cutoff_score
    tied <- ordered_scores == cutoff_score
    n_above <- sum(above)
    slots_in_tie <- k - n_above
    expected_tp <- sum(ordered_labels[above]) +
      slots_in_tie * mean(ordered_labels[tied])
    prec <- expected_tp / k
    rec <- expected_tp / n_positive
    data.frame(
      top_k_requested = k,
      n_selected = k,
      cutoff_score = cutoff_score,
      precision = prec,
      recall = rec,
      f1 = if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0
    )
  })

  list(
    auroc = auroc,
    average_precision = average_precision,
    trapezoidal_auprc = trapezoidal_auprc,
    topk = do.call(rbind, topk),
    order = ord,
    block_start = block_start,
    block_end = block_end
  )
}

threshold_average_precision <- function(labels, scores) {
  labels <- as.logical(labels)
  keep <- !is.na(labels) & is.finite(scores)
  labels <- labels[keep]
  scores <- scores[keep]
  positives <- sum(labels)
  if (!length(labels) || positives == 0L) {
    return(NA_real_)
  }
  ordering <- order(scores, decreasing = TRUE)
  labels <- labels[ordering]
  scores <- scores[ordering]
  blocks <- rle(scores)$lengths
  block_end <- cumsum(blocks)
  block_start <- c(1L, head(block_end, -1L) + 1L)
  block_positives <- vapply(seq_along(blocks), function(index) {
    sum(labels[block_start[[index]]:block_end[[index]]])
  }, numeric(1))
  precision <- cumsum(block_positives) / block_end
  sum((block_positives / positives) * precision)
}

tie_safe_auroc <- function(labels, scores) {
  labels <- as.logical(labels)
  keep <- !is.na(labels) & is.finite(scores)
  labels <- labels[keep]
  scores <- scores[keep]
  positives <- sum(labels)
  negatives <- length(labels) - positives
  if (positives == 0L || negatives == 0L) {
    return(NA_real_)
  }
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[labels]) - positives * (positives + 1) / 2) /
    (as.double(positives) * negatives)
}

expected_precision_at_k <- function(labels, scores, k) {
  labels <- as.logical(labels)
  keep <- !is.na(labels) & is.finite(scores)
  labels <- labels[keep]
  scores <- scores[keep]
  k <- min(as.integer(k), length(labels))
  if (!length(labels) || k <= 0L) {
    return(NA_real_)
  }
  blocks <- data.table::data.table(label = as.integer(labels), score = scores)[
    , .(size = .N, positives = sum(label)),
    by = score
  ][order(-score)]
  remaining <- k
  expected_true <- 0
  for (index in seq_len(nrow(blocks))) {
    if (remaining <= 0L) break
    selected <- min(remaining, blocks$size[[index]])
    expected_true <- expected_true +
      selected * blocks$positives[[index]] / blocks$size[[index]]
    remaining <- remaining - selected
  }
  expected_true / k
}

finite_mean <- function(values) {
  keep <- is.finite(values)
  if (!any(keep)) {
    return(NA_real_)
  }
  mean(values[keep])
}

perturbation_response_stat <- function(x) {
  p <- pmax(as.numeric(x$PValue), .Machine$double.xmin)
  stat <- abs(as.numeric(x$logFC)) * -log10(p)
  stat[!is.finite(stat)] <- max(stat[is.finite(stat)], 0) + 1
  names(stat) <- x$feature
  stat <- stat[!duplicated(names(stat))]
  stat <- sort(stat + seq_along(stat) * 1e-14, decreasing = TRUE)
  stat
}
