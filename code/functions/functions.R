prepare_ground_truth <- function(ground_truth) {
  if (is.null(ground_truth) || ncol(ground_truth) < 2) {
    stop("ground_truth must contain at least two columns.")
  }

  gt <- ground_truth
  gt_names <- tolower(colnames(gt))

  if (all(c("regulator", "target") %in% gt_names)) {
    regulator_col <- which(gt_names == "regulator")[1]
    target_col <- which(gt_names == "target")[1]
  } else if (all(c("gene1", "gene2") %in% gt_names)) {
    regulator_col <- which(gt_names == "gene1")[1]
    target_col <- which(gt_names == "gene2")[1]
  } else if (identical(gt_names[1:2], c("gene2", "gene1"))) {
    regulator_col <- 2
    target_col <- 1
  } else {
    regulator_col <- 1
    target_col <- 2
  }

  gt <- gt[, c(regulator_col, target_col), drop = FALSE]
  colnames(gt) <- c("regulator", "target")
  gt
}

prepare_network_table <- function(network_table) {
  if (is.null(network_table) || ncol(network_table) < 2) {
    stop("network_table must contain at least two columns.")
  }

  nt <- as.data.frame(network_table, stringsAsFactors = FALSE)
  nt_names <- tolower(colnames(nt))

  if (all(c("regulator", "target") %in% nt_names)) {
    regulator_col <- which(nt_names == "regulator")[1]
    target_col <- which(nt_names == "target")[1]
  } else if (all(c("gene1", "gene2") %in% nt_names)) {
    regulator_col <- which(nt_names == "gene1")[1]
    target_col <- which(nt_names == "gene2")[1]
  } else {
    regulator_col <- 1
    target_col <- 2
  }

  weight_col <- NA_integer_
  weight_candidates <- c(
    "weight",
    "score",
    "cor_value",
    "importance",
    "corr",
    "value"
  )
  for (candidate in weight_candidates) {
    if (candidate %in% nt_names) {
      weight_col <- which(nt_names == candidate)[1]
      break
    }
  }
  if (is.na(weight_col)) {
    if (ncol(nt) >= 3) {
      weight_col <- 3
    } else {
      nt$weight <- 1
      weight_col <- ncol(nt)
    }
  }

  nt <- nt[, c(regulator_col, target_col, weight_col), drop = FALSE]
  colnames(nt) <- c("regulator", "target", "weight")
  nt$regulator <- as.character(nt$regulator)
  nt$target <- as.character(nt$target)
  nt$weight <- suppressWarnings(as.numeric(nt$weight))
  nt <- nt[complete.cases(nt[, c("regulator", "target")]), , drop = FALSE]
  nt
}

square_matrix_to_network <- function(
  weight_matrix,
  remove_self = TRUE
) {
  mat <- as.matrix(weight_matrix)

  if (is.null(rownames(mat))) {
    rownames(mat) <- seq_len(nrow(mat))
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- seq_len(ncol(mat))
  }

  index <- which(!is.na(mat), arr.ind = TRUE)
  if (nrow(index) == 0) {
    return(data.frame(
      regulator = character(0),
      target = character(0),
      weight = numeric(0)
    ))
  }

  network_table <- data.frame(
    regulator = colnames(mat)[index[, 2]],
    target = rownames(mat)[index[, 1]],
    weight = as.numeric(mat[index]),
    stringsAsFactors = FALSE
  )

  if (remove_self) {
    network_table <- network_table[
      network_table$regulator != network_table$target, ,
      drop = FALSE
    ]
  }

  network_table <- network_table[
    order(
      abs(network_table$weight),
      decreasing = TRUE
    ), ,
    drop = FALSE
  ]
  rownames(network_table) <- NULL
  network_table
}

safe_metric_divide <- function(numerator, denominator, default = 0) {
  if (is.na(denominator) || denominator == 0) {
    return(default)
  }
  numerator / denominator
}

# BEELINE EP: select all nonzero predictions tied at the top-k boundary.
# Unlike a fixed edge-budget statistic, divide by the actual selected count.
select_early_precision_edges <- function(weights, top_k) {
  eligible <- which(is.finite(weights) & weights > 0)
  if (!length(eligible) || top_k <= 0) return(integer(0))
  ranked <- eligible[order(-weights[eligible])]
  cutoff <- weights[ranked[[min(as.integer(top_k), length(ranked))]]]
  eligible[weights[eligible] >= cutoff]
}

metric_edge_ids <- function(regulator, target) {
  regulator <- enc2utf8(as.character(regulator))
  target <- enc2utf8(as.character(target))
  paste0(nchar(regulator, type = "bytes"), ":", regulator, target, recycle0 = TRUE)
}

fast_standard_metric_summary <- function(network_table, ground_truth) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(
      "Package 'data.table' is required for fast standard metric evaluation."
    )
  }

  gt <- prepare_ground_truth(ground_truth)
  gt <- gt[gt$regulator != gt$target, , drop = FALSE]
  gt <- unique(gt)
  if (!nrow(gt)) {
    return(list(
      AUROC = NA_real_,
      AUPRC = NA_real_,
      ACC = NA_real_,
      Precision = NA_real_,
      Recall = NA_real_,
      F1 = NA_real_,
      JI = NA_real_,
      SI = 0,
      EarlyPrecision = NA_real_,
      EPR = NA_real_,
      TopK = 0L,
      EarlySelectedEdges = 0L,
      EarlyTruePositives = 0L,
      PositiveRate = NA_real_,
      GoldEdges = 0L
    ))
  }

  pred <- prepare_network_table(network_table)
  gt_genes <- sort(unique(c(gt$regulator, gt$target)))
  pred <- pred[
    pred$regulator %in%
      gt_genes &
      pred$target %in% gt_genes &
      pred$regulator != pred$target, ,
    drop = FALSE
  ]

  gt_dt <- data.table::as.data.table(gt)
  gt_dt[, edge_id := metric_edge_ids(regulator, target)]
  gt_dt <- unique(gt_dt, by = "edge_id")
  true_edge_ids <- gt_dt$edge_id

  pred_dt <- data.table::as.data.table(pred)
  pred_dt[, regulator := as.character(regulator)]
  pred_dt[, target := as.character(target)]
  pred_dt[, weight := suppressWarnings(as.numeric(weight))]
  pred_dt <- pred_dt[is.finite(weight)]
  pred_dt[, weight := abs(weight)]
  pred_dt[, edge_id := metric_edge_ids(regulator, target)]
  data.table::setorder(pred_dt, -weight)
  pred_dt <- pred_dt[!duplicated(edge_id)]
  pred_dt[, label := as.integer(edge_id %in% true_edge_ids)]

  n_genes <- length(gt_genes)
  n_true <- length(true_edge_ids)
  n_possible <- as.double(n_genes) * as.double(max(n_genes - 1L, 0L))
  n_negative <- n_possible - n_true

  explicit_total <- nrow(pred_dt)
  explicit_pos <- sum(pred_dt$label)
  explicit_neg <- explicit_total - explicit_pos
  missing_pos <- n_true - explicit_pos
  missing_neg <- n_negative - explicit_neg
  if (missing_pos < 0 || missing_neg < 0) {
    stop("Fast metric summary encountered an invalid edge accounting state.")
  }

  grouped <- pred_dt[,
    .(
      pos = sum(label == 1L),
      neg = sum(label == 0L)
    ),
    by = .(weight)
  ]

  if (missing_pos > 0 || missing_neg > 0) {
    zero_idx <- which(grouped$weight == 0)
    if (length(zero_idx)) {
      grouped$pos[zero_idx[[1]]] <- grouped$pos[zero_idx[[1]]] + missing_pos
      grouped$neg[zero_idx[[1]]] <- grouped$neg[zero_idx[[1]]] + missing_neg
    } else {
      grouped <- data.table::rbindlist(
        list(
          grouped,
          data.table::data.table(
            weight = 0,
            pos = missing_pos,
            neg = missing_neg
          )
        ),
        use.names = TRUE
      )
    }
  }

  data.table::setorder(grouped, -weight)
  tp_cum <- cumsum(grouped$pos)
  fp_cum <- cumsum(grouped$neg)

  positive_rate <- safe_metric_divide(n_true, n_possible, default = NA_real_)
  top_k <- as.integer(n_true)
  selected <- select_early_precision_edges(pred_dt$weight, top_k)
  top_n <- length(selected)
  top_hits <- sum(pred_dt$label[selected])
  early_precision_raw <- safe_metric_divide(
    top_hits,
    top_n,
    default = 0
  )
  epr_raw <- safe_metric_divide(
    early_precision_raw,
    positive_rate,
    default = NA_real_
  )

  if (n_true == 0 || n_negative <= 0) {
    auroc_raw <- NA_real_
  } else {
    grouped_asc <- grouped[order(weight)]
    neg_before <- c(0, head(cumsum(grouped_asc$neg), -1))
    auroc_raw <- safe_metric_divide(
      sum(grouped_asc$pos * (neg_before + 0.5 * grouped_asc$neg)),
      n_true * n_negative,
      default = NA_real_
    )
  }

  if (n_true == 0) {
    auprc_raw <- NA_real_
  } else {
    precision_desc <- ifelse(
      (tp_cum + fp_cum) > 0,
      tp_cum / (tp_cum + fp_cum),
      0
    )
    recall_desc <- tp_cum / n_true
    rev_idx <- rev(seq_along(precision_desc))
    pr_precision <- c(precision_desc[rev_idx], 1)
    pr_recall <- c(recall_desc[rev_idx], 0)
    auprc_raw <- abs(sum(
      diff(pr_recall) * (head(pr_precision, -1) + tail(pr_precision, -1)) / 2
    ))
  }

  if (n_negative <= 0) {
    tp_best <- 0
    fp_best <- 0
  } else {
    # Compare count numerators of J=TP/P-FP/N. which.max takes the first
    # (highest-score) threshold when the mathematical maxima are tied.
    best_idx <- which.max(as.double(tp_cum) * n_negative - as.double(fp_cum) * n_true)
    tp_best <- tp_cum[[best_idx]]
    fp_best <- fp_cum[[best_idx]]
  }
  tn_best <- n_negative - fp_best
  fn_best <- n_true - tp_best
  precision_best <- safe_metric_divide(tp_best, tp_best + fp_best, default = 0)
  recall_best <- safe_metric_divide(tp_best, tp_best + fn_best, default = 0)
  f1_best <- if ((precision_best + recall_best) == 0) {
    0
  } else {
    2 * precision_best * recall_best / (precision_best + recall_best)
  }
  acc_best <- safe_metric_divide(
    tp_best + tn_best,
    n_true + n_negative,
    default = 0
  )
  union_n <- (tp_best + fp_best) + n_true - tp_best
  ji_best <- if (union_n == 0) 0 else tp_best / union_n

  list(
    AUROC = as.numeric(auroc_raw),
    AUPRC = as.numeric(auprc_raw),
    ACC = as.numeric(acc_best),
    Precision = as.numeric(precision_best),
    Recall = as.numeric(recall_best),
    F1 = as.numeric(f1_best),
    JI = as.numeric(ji_best),
    SI = as.integer(tp_best),
    EarlyPrecision = early_precision_raw,
    EPR = as.numeric(epr_raw),
    TopK = top_k,
    EarlySelectedEdges = top_n,
    EarlyTruePositives = as.integer(top_hits),
    PositiveRate = positive_rate,
    GoldEdges = as.integer(n_true)
  )
}
