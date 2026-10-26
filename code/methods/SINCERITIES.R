suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(ppcor))

# Method parameters are passed as function arguments.

sincerities_use_upstream <- function(use_upstream = FALSE) {
  isTRUE(use_upstream)
}

sincerities_parallel_verbose <- function(verbose = FALSE) {
  isTRUE(verbose)
}

sincerities_max_time_points <- function(max_time_points = NA_integer_) {
  if (is.null(max_time_points) || is.na(max_time_points)) {
    return(NA_integer_)
  }
  parsed <- suppressWarnings(as.integer(max_time_points))
  if (!is.finite(parsed) || parsed <= 1L) {
    return(NA_integer_)
  }
  as.integer(parsed)
}

sincerities_ranked_predictions <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- dget(file.path(benchmark_root_dir(), "vendor", "SINCERITIES", "final_ranked_predictions.R"))
    }
    cached
  }
})

sincerities_parallel_map <- function(x, fun, cores = 1L, export_fun = character(0), verbose = FALSE) {
  if (cores <= 1L || !requireNamespace("thisutils", quietly = TRUE)) {
    return(lapply(x, fun))
  }

  thisutils::parallelize_fun(
    x = x,
    fun = fun,
    cores = cores,
    export_fun = export_fun,
    throw_error = TRUE,
    clean_result = FALSE,
    verbose = verbose
  )
}

compress_sincerities_pseudotime <- function(pseudotime, max_time_points = sincerities_max_time_points()) {
  pseudotime <- as.numeric(pseudotime)
  if (!is.finite(max_time_points) || is.na(max_time_points) || max_time_points <= 1L) {
    return(pseudotime)
  }

  unique_count <- length(unique(pseudotime))
  if (unique_count <= max_time_points) {
    return(pseudotime)
  }

  ranks <- rank(pseudotime, ties.method = "first", na.last = "keep")
  keep <- is.finite(ranks)
  binned <- pseudotime
  binned[keep] <- floor((ranks[keep] - 1) * max_time_points / sum(keep))
  as.numeric(binned)
}

build_sincerities_input <- function(matrix, pseudotime, max_time_points = sincerities_max_time_points()) {
  expression_matrix <- t(as.matrix(matrix))

  if (is.data.frame(pseudotime) || is.matrix(pseudotime)) {
    pseudotime <- pseudotime[, 1]
  }
  pseudotime <- as.numeric(pseudotime)

  if (length(pseudotime) != nrow(expression_matrix)) {
    stop("Length of pseudotime does not match the number of cells.")
  }

  if (!is.null(names(pseudotime)) && all(rownames(expression_matrix) %in% names(pseudotime))) {
    pseudotime <- pseudotime[rownames(expression_matrix)]
  }

  pseudotime <- compress_sincerities_pseudotime(
    pseudotime = pseudotime,
    max_time_points = max_time_points
  )

  ordered_time <- sort(unique(pseudotime))
  single_cell_data <- lapply(ordered_time, function(timepoint) {
    expression_matrix[pseudotime == timepoint, , drop = FALSE]
  })

  list(
    time = ordered_time,
    num_time_points = length(ordered_time),
    totDATA = expression_matrix,
    timeline = pseudotime,
    numGENES = ncol(expression_matrix),
    genes = colnames(expression_matrix),
    singleCELLdata = single_cell_data
  )
}

sincerities_nonconstant_gene_mask <- function(matrix) {
  if (is.null(matrix) || !nrow(matrix)) {
    return(logical(0))
  }
  apply(matrix, 1, function(values) {
    values <- suppressWarnings(as.numeric(values))
    keep <- is.finite(values)
    if (sum(keep) < 2) {
      return(FALSE)
    }
    stats::sd(values[keep]) > 0
  })
}

run_SINCERITIES_parallel <- function(
  DATA,
  method = 1,
  noDIAG = 1,
  SIGN = 0,
  cores = 1L,
  verbose = FALSE
) {
  if (!method %in% c(1, 2, 3, 4)) {
    stop("Choose regularization regression strategy with 1,2,3,4.")
  }
  if (!noDIAG %in% c(0, 1)) {
    stop("noDIAG should be either 0 or 1")
  }
  if (!SIGN %in% c(0, 1)) {
    stop("SIGN should be either 0 or 1")
  }

  single_cell_data <- DATA$singleCELLdata
  time <- DATA$time
  numGENES <- DATA$numGENES
  num_time_points <- length(time)

  if (num_time_points < 5) {
    stop("** DATA with a number of time points < 5. Please run SINCERITIES_CROSS_VALIDATION function **")
  }

  distance_cols <- sincerities_parallel_map(
    x = seq_len(numGENES),
    cores = cores,
    verbose = verbose,
    fun = function(gi) {
      values <- numeric(num_time_points - 1)
      for (ti in seq_len(num_time_points - 1)) {
        data_ti <- t(single_cell_data[[ti]])
        data_ti_plus1 <- t(single_cell_data[[ti + 1]])
        p1 <- data_ti[gi, ]
        p2 <- data_ti_plus1[gi, ]
        values[ti] <- ks.test(p1, p2)$statistic
      }
      values
    }
  )

  DISTANCE_matrix <- do.call(cbind, distance_cols)
  deltaT <- replicate(ncol(DISTANCE_matrix), time[2:length(time)] - time[1:(length(time) - 1)])
  DISTANCE_matrix <- DISTANCE_matrix / deltaT

  if (method == 1) {
    alphas <- 0
  } else if (method == 2) {
    alphas <- seq(0, 1, 0.1)
  } else if (method == 3) {
    alphas <- 1
  } else {
    stop("Manual alpha selection is not supported in benchmark mode.")
  }

  X_matrix <- DISTANCE_matrix[1:(num_time_points - 2), , drop = FALSE]
  predictor_keep_global <- apply(X_matrix, 2, function(values) {
    values <- suppressWarnings(as.numeric(values))
    keep <- is.finite(values)
    if (sum(keep) < 2) {
      return(FALSE)
    }
    stats::sd(values[keep]) > 0
  })
  nfold <- nrow(X_matrix)
  foldid <- seq_len(nfold)

  gene_results <- sincerities_parallel_map(
    x = seq_len(numGENES),
    cores = cores,
    verbose = verbose,
    fun = function(gi) {
      suppressPackageStartupMessages(library(glmnet))
      Y_vector <- DISTANCE_matrix[2:(num_time_points - 1), gi]
      full_beta <- numeric(ncol(X_matrix))
      y_keep <- is.finite(Y_vector)
      if (sum(y_keep) < 2 || stats::sd(Y_vector[y_keep]) == 0) {
        return(list(
          lambda = NA_real_,
          alpha = NA_real_,
          beta = full_beta
        ))
      }

      predictor_keep <- predictor_keep_global
      if (noDIAG == 1 && gi <= length(predictor_keep)) {
        predictor_keep[gi] <- FALSE
      }
      valid_cols <- which(predictor_keep)
      if (!length(valid_cols)) {
        return(list(
          lambda = NA_real_,
          alpha = NA_real_,
          beta = full_beta
        ))
      }

      X_fit <- X_matrix[, valid_cols, drop = FALSE]
      lambda <- numeric(length(alphas))
      cvERROR <- numeric(length(alphas))
      beta <- matrix(data = 0, nrow = ncol(X_matrix), ncol = length(alphas))

      for (test in seq_along(alphas)) {
        CV_results <- tryCatch(
          cv.glmnet(
            X_fit,
            Y_vector,
            alpha = alphas[test],
            nfolds = nfold,
            foldid = foldid,
            keep = TRUE,
            lower.limits = 0,
            upper.limits = Inf,
            grouped = FALSE
          ),
          error = function(...) NULL
        )
        if (is.null(CV_results)) {
          lambda[test] <- NA_real_
          cvERROR[test] <- Inf
          next
        }
        lambda[test] <- CV_results$lambda.min
        cvERROR[test] <- CV_results$cvm[CV_results$lambda == CV_results$lambda.min][1]
        coef.CV_results <- coef(CV_results, s = "lambda.min")
        if (length(coef.CV_results@i) > 1) {
          selected_idx <- coef.CV_results@i[-1]
          beta[valid_cols[selected_idx], test] <- coef.CV_results@x[-1]
        }
      }

      if (!any(is.finite(cvERROR))) {
        return(list(
          lambda = NA_real_,
          alpha = NA_real_,
          beta = full_beta
        ))
      }

      minIdx <- max(which(cvERROR == min(cvERROR)))
      list(
        lambda = lambda[minIdx],
        alpha = alphas[minIdx],
        beta = beta[, minIdx]
      )
    }
  )

  pred_lambda_min <- do.call(cbind, lapply(gene_results, `[[`, "beta"))

  if (SIGN == 1) {
    parcorr_matrix <- pcor(DATA$totDATA, method = "spearman")$estimate
    pred_lambda_min <- pred_lambda_min * sign(parcorr_matrix)
  }

  list(
    DISTANCE_matrix = DISTANCE_matrix,
    adj_matrix = pred_lambda_min,
    lambda = vapply(gene_results, `[[`, numeric(1), "lambda"),
    alpha = vapply(gene_results, `[[`, numeric(1), "alpha")
  )
}

run_SINCERITIES <- function(
  matrix,
  pseudotime,
  regulators = NULL,
  targets = NULL,
  filename = "",
  cores = 1L,
  use_upstream = sincerities_use_upstream(),
  verbose = sincerities_parallel_verbose()
) {
  gene_keep <- rownames(matrix)
  if (!is.null(regulators) && length(regulators) > 0) {
    gene_keep <- intersect(gene_keep, regulators)
  }
  if (!is.null(targets) && length(targets) > 0) {
    gene_keep <- union(gene_keep, intersect(rownames(matrix), targets))
  }
  matrix <- matrix[unique(gene_keep), , drop = FALSE]
  nonconstant_mask <- sincerities_nonconstant_gene_mask(matrix)
  matrix <- matrix[nonconstant_mask, , drop = FALSE]
  if (nrow(matrix) < 2) {
    return(NULL)
  }

  network_table <- try({
    expression_data <- build_sincerities_input(matrix, pseudotime)
    result <- run_SINCERITIES_parallel(
      DATA = expression_data,
      method = 1,
      noDIAG = 1,
      SIGN = 0,
      cores = if (isTRUE(use_upstream)) 1L else cores,
      verbose = verbose
    )
    adj_matrix <- result$adj_matrix
    max_value <- suppressWarnings(max(abs(adj_matrix), na.rm = TRUE))
    if (is.finite(max_value) && max_value > 0) {
      adj_matrix <- adj_matrix / max_value
    }
    network_table <- as.data.frame(
      sincerities_ranked_predictions()(
        adj_matrix,
        expression_data$genes,
        SIGN = 1,
        saveFile = FALSE
      )[, 1:3]
    )
    colnames(network_table) <- c("regulator", "target", "weight")
    if (!is.null(regulators) && length(regulators) > 0) {
      network_table <- network_table[network_table$regulator %in% regulators, , drop = FALSE]
    }
    if (!is.null(targets) && length(targets) > 0) {
      network_table <- network_table[network_table$target %in% targets, , drop = FALSE]
    }
    network_table
  })

  if (inherits(network_table, "try-error")) {
    return(NULL)
  }

  write.csv(
    network_table,
    file = filename,
    row.names = FALSE
  )

  network_table
}
