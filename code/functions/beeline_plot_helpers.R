has_complex_heatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
  requireNamespace("circlize", quietly = TRUE) &&
  requireNamespace("grid", quietly = TRUE)

heatmap_palette <- function(metric_name) {
  palette_name <- switch(as.character(metric_name),
    AUPRC = "viridis",
    EarlyPrecision = "plasma",
    EPR = "Spectral",
    AUROC = "RdBu",
    ACC = "magma",
    Precision = "inferno",
    Recall = "cividis",
    F1 = "RdYlBu",
    "RdBu"
  )

  if (
    palette_name %in%
      c("viridis", "plasma", "magma", "inferno", "cividis") &&
      requireNamespace("viridisLite", quietly = TRUE)
  ) {
    return(viridisLite::viridis(11, option = palette_name, direction = 1))
  }
  if (
    requireNamespace("RColorBrewer", quietly = TRUE) &&
      palette_name %in% rownames(RColorBrewer::brewer.pal.info)
  ) {
    max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
    return(grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(
      max_colors,
      palette_name
    )))(11))
  }
  grDevices::hcl.colors(11, "Blue-Red 3")
}

heatmap_colors <- function(values, metric_name = NULL) {
  rng <- range(values, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    rng <- c(0, 1)
  }
  colors <- heatmap_palette(metric_name)
  circlize::colorRamp2(
    seq(rng[[1]], rng[[2]], length.out = length(colors)),
    colors
  )
}

heatmap_legend_style <- function() {
  list(
    title_gp = grid::gpar(fontsize = 7, fontface = "plain"),
    labels_gp = grid::gpar(fontsize = 7, fontface = "plain"),
    grid_width = grid::unit(0.118, "in")
  )
}

significance_label_beeline <- function(p_value) {
  if (is.na(p_value) || p_value > 0.05) {
    return(NA_character_)
  }
  if (p_value > 0.01) {
    return("*")
  }
  if (p_value > 0.001) {
    return("**")
  }
  if (p_value > 1e-04) {
    return("***")
  }
  "****"
}

normalize_beeline_method_name <- function(method_name) {
  normalized <- as.character(method_name)
  normalized[normalized %in% c("SINCERITITES", "SINCERITIES")] <- "SINCERITIES"
  normalized
}

preferred_beeline_method_order <- function(methods) {
  methods <- unique(as.character(methods))
  front <- intersect("multiCSN", methods)
  c(front, sort(setdiff(methods, front)))
}

beeline_method_colors <- function(methods) {
  fixed_colors <- c(
    "multiCSN" = "#D70440",
    "inferCSN" = "#6A3D9A",
    "abess" = "#4575B4",
    "abess_ebic" = "#74ADD1",
    "GENIE3" = "#2177B8",
    "GRNBoost2" = "#5E4FA2",
    "PIDC" = "#2C7FB8",
    "PPCOR" = "#0AA344",
    "LEAP" = "#F9BD10",
    "SINCERITIES" = "#8076A3",
    "SCODE" = "#CC6F00"
  )

  methods <- unique(methods)
  missing_methods <- setdiff(methods, names(fixed_colors))
  if (length(missing_methods)) {
    extra_colors <- grDevices::hcl.colors(
      n = length(missing_methods),
      palette = "Dynamic"
    )
    names(extra_colors) <- missing_methods
    fixed_colors <- c(fixed_colors, extra_colors)
  }

  fixed_colors[methods]
}

beeline_method_shapes <- function(methods) {
  fixed_shapes <- c(16, 17, 15, 18, 3, 7, 8, 1, 2, 0, 4, 9)
  methods <- unique(methods)
  shape_values <- rep(fixed_shapes, length.out = length(methods))
  names(shape_values) <- methods
  shape_values
}

parse_synthetic_dataset_info <- function(dataset_values) {
  split_values <- strsplit(as.character(dataset_values), "-", fixed = TRUE)
  family <- vapply(
    split_values,
    function(parts) {
      if (length(parts) <= 2) {
        return(parts[[1]])
      }
      paste(parts[seq_len(length(parts) - 2)], collapse = "-")
    },
    character(1)
  )
  cells <- vapply(
    split_values,
    function(parts) {
      if (length(parts) >= 2) {
        parts[[length(parts) - 1]]
      } else {
        NA_character_
      }
    },
    character(1)
  )
  replicate_id <- vapply(
    split_values,
    function(parts) {
      if (length(parts) >= 1) {
        parts[[length(parts)]]
      } else {
        NA_character_
      }
    },
    character(1)
  )

  data.frame(
    Family = family,
    Cells = cells,
    Replicate = replicate_id,
    stringsAsFactors = FALSE
  )
}

metric_is_higher_better <- function(metric_name) {
  !identical(metric_name, "RuntimeSeconds")
}

build_matrix_from_df <- function(df, row_col, col_col, value_col) {
  rows <- unique(as.character(df[[row_col]]))
  cols <- unique(as.character(df[[col_col]]))
  mat <- matrix(
    NA_real_,
    nrow = length(rows),
    ncol = length(cols),
    dimnames = list(rows, cols)
  )
  for (i in seq_len(nrow(df))) {
    r <- as.character(df[[row_col]][[i]])
    c <- as.character(df[[col_col]][[i]])
    v <- suppressWarnings(as.numeric(df[[value_col]][[i]]))
    if (!is.finite(v)) {
      next
    }
    mat[r, c] <- v
  }
  mat
}

multicsn_signif_labels_from_matrix <- function(
  mat,
  metric_name = NULL,
  method_order = colnames(mat),
  baseline_method = "multiCSN"
) {
  if (is.null(mat) || !nrow(mat) || !ncol(mat)) {
    return(character(0))
  }
  if (is.null(method_order) || !length(method_order)) {
    method_order <- colnames(mat)
  }
  labels <- rep("", length(method_order))
  names(labels) <- method_order
  if (is.null(colnames(mat)) || !(baseline_method %in% colnames(mat))) {
    return(labels)
  }

  alt <- if (
    !is.null(metric_name) &&
      nzchar(metric_name) &&
      !metric_is_higher_better(metric_name)
  ) {
    "less"
  } else {
    "greater"
  }

  baseline_values <- suppressWarnings(as.numeric(mat[, baseline_method]))
  baseline_values <- baseline_values[is.finite(baseline_values)]
  if (!length(baseline_values)) {
    return(labels)
  }

  for (method in method_order) {
    if (identical(method, baseline_method)) {
      next
    }
    if (!method %in% colnames(mat)) {
      next
    }
    other_values <- suppressWarnings(as.numeric(mat[, method]))
    other_values <- other_values[is.finite(other_values)]
    if (!length(other_values)) {
      next
    }
    p_value <- tryCatch(
      stats::wilcox.test(
        baseline_values,
        other_values,
        alternative = alt,
        exact = FALSE
      )$p.value,
      error = function(...) NA_real_
    )
    lab <- significance_label_beeline(p_value)
    if (!is.na(lab)) {
      labels[[method]] <- lab
    }
  }
  labels
}

make_metric_heatmap <- function(
  mat,
  title,
  metric_name = NULL,
  row_split = NULL,
  col_order = NULL,
  legend_title = NULL,
  method_meta = NULL,
  cell_size_in = 0.157,
  show_box_axis = TRUE,
  show_legend = TRUE,
  show_significance = TRUE
) {
  if (
    !isTRUE(has_complex_heatmap) || is.null(mat) || !nrow(mat) || !ncol(mat)
  ) {
    return(NULL)
  }

  if (!is.null(col_order)) {
    keep <- intersect(col_order, colnames(mat))
    mat <- mat[, keep, drop = FALSE]
  }

  method_order <- colnames(mat)
  method_colors <- NULL
  if (!is.null(method_meta) && !is.null(method_meta$colors)) {
    method_colors <- unname(method_meta$colors[method_order])
  }
  if (is.null(method_colors) || any(is.na(method_colors))) {
    method_colors <- rep("#CBD5E1", length(method_order))
  }

  sig_labels <- if (isTRUE(show_significance)) {
    multicsn_signif_labels_from_matrix(
      mat = mat,
      metric_name = metric_name,
      method_order = method_order,
      baseline_method = "multiCSN"
    )
  } else {
    stats::setNames(rep("", length(method_order)), method_order)
  }
  sig_labels <- sig_labels[method_order]
  sig_labels[is.na(sig_labels)] <- ""

  box_values <- lapply(method_order, function(method) {
    vals <- suppressWarnings(as.numeric(mat[, method]))
    vals <- vals[is.finite(vals)]
    if (!length(vals)) {
      return(0)
    }
    vals
  })
  names(box_values) <- method_order

  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    sig = ComplexHeatmap::anno_text(
      sig_labels,
      gp = grid::gpar(col = "gray30", fontsize = 10),
      just = "center"
    ),
    Distribution = ComplexHeatmap::anno_boxplot(
      box_values,
      gp = grid::gpar(fill = method_colors, col = "gray30", lwd = 0.5),
      outline = TRUE,
      axis = isTRUE(show_box_axis),
      border = FALSE
    ),
    show_annotation_name = FALSE,
    annotation_height = grid::unit.c(
      grid::unit(0.138, "in"),
      grid::unit(0.866, "in")
    )
  )

  col_fun <- heatmap_colors(as.numeric(mat), metric_name)
  ComplexHeatmap::Heatmap(
    mat,
    name = if (!is.null(legend_title)) legend_title else "Value",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = row_split,
    column_title = title,
    row_names_side = "left",
    show_row_names = TRUE,
    top_annotation = top_anno,
    show_heatmap_legend = isTRUE(show_legend),
    heatmap_legend_param = heatmap_legend_style(),
    width = grid::unit(ncol(mat) * cell_size_in, "in"),
    height = grid::unit(nrow(mat) * cell_size_in, "in")
  )
}

plot_heatmaps_beeline <- function(
  evaluation_data,
  dataset_type
) {
  heatmaps <- list()
  metrics <- resolve_plot_metrics(evaluation_data)
  if (!length(metrics) || !isTRUE(has_complex_heatmap)) {
    return(invisible(NULL))
  }
  method_meta <- resolve_plot_methods(evaluation_data)
  df <- evaluation_data
  df <- df[df$Method %in% method_meta$order, , drop = FALSE]
  if (!nrow(df)) {
    return(invisible(NULL))
  }
  df$Method <- factor(df$Method, levels = method_meta$order)

  if (identical(dataset_type, "synthetic")) {
    cell_med <- synthetic_cell_level_medians(df, metrics)
    if (is.null(cell_med) || !nrow(cell_med)) {
      return(invisible(NULL))
    }
    families <- levels(cell_med$Family)

    cell_med$RowKey <- paste0(
      as.character(cell_med$Family),
      "-",
      as.character(cell_med$Cells)
    )
    row_order <- unique(cell_med$RowKey[order(
      as.character(cell_med$Family),
      cell_med$CellsNum
    )])
    row_split_map <- setNames(
      as.character(cell_med$Family),
      as.character(cell_med$RowKey)
    )
    row_split <- factor(row_split_map[row_order], levels = families)

    mats_sub <- list()
    for (metric in metrics) {
      sub <- cell_med[cell_med$Metric == metric, , drop = FALSE]
      if (!nrow(sub)) {
        next
      }
      mat_df <- sub |>
        dplyr::group_by(RowKey, Method) |>
        dplyr::summarise(
          Value = stats::median(Value, na.rm = TRUE),
          .groups = "drop"
        )
      mat <- build_matrix_from_df(mat_df, "RowKey", "Method", "Value")
      mat <- mat[row_order, method_meta$order, drop = FALSE]
      mats_sub[[metric]] <- mat
    }
    if (!length(mats_sub)) {
      return(invisible(NULL))
    }

    cell_size_in <- 0.157
    gap_in <- 0.276

    ht_sub <- NULL
    for (i in seq_along(metrics)) {
      metric <- metrics[[i]]
      ht <- make_metric_heatmap(
        mat = mats_sub[[metric]],
        title = metric,
        metric_name = metric,
        row_split = row_split,
        col_order = method_meta$order,
        legend_title = metric,
        method_meta = method_meta,
        cell_size_in = cell_size_in,
        show_box_axis = TRUE,
        show_legend = TRUE,
        show_significance = FALSE
      )
      if (is.null(ht)) {
        next
      }
      ht_sub <- if (is.null(ht_sub)) ht else (ht_sub + ht)
    }
    if (!is.null(ht_sub)) {
      heatmaps[["subdataset"]] <- list(
        heatmap = ht_sub,
        gap_in = gap_in
      )
    }


  }

  if (identical(dataset_type, "curated")) {
    dataset_parts <- strsplit(as.character(df$Dataset), "-", fixed = TRUE)
    drop_pct <- vapply(
      seq_along(dataset_parts),
      function(i) {
        parts <- dataset_parts[[i]]
        if (length(parts) >= 4) {
          as.numeric(parts[[length(parts)]])
        } else {
          0
        }
      },
      numeric(1)
    )
    drop_pct[!is.finite(drop_pct)] <- 0
    df$DropoutPct <- drop_pct
    df$DropoutSuffix <- ifelse(
      df$DropoutPct > 0,
      paste0("-", df$DropoutPct),
      ""
    )

    cell_size_in <- 0.157
    gap_in <- 0.276

    mats_dropout <- list()
    dropout_levels <- sort(unique(df$DropoutPct))
    dropout_labels <- as.character(dropout_levels)
    for (metric in metrics) {
      mat_df <- df |>
        dplyr::group_by(DropoutPct, Method) |>
        dplyr::summarise(
          Value = stats::median(.data[[metric]], na.rm = TRUE),
          .groups = "drop"
        )
      mat_df$RowKey <- as.character(mat_df$DropoutPct)
      mat <- build_matrix_from_df(mat_df, "RowKey", "Method", "Value")
      keep_rows <- intersect(dropout_labels, rownames(mat))
      mat <- mat[keep_rows, method_meta$order, drop = FALSE]
      mats_dropout[[metric]] <- mat
    }
    ht_dropout <- NULL
    for (i in seq_along(metrics)) {
      metric <- metrics[[i]]
      ht <- make_metric_heatmap(
        mat = mats_dropout[[metric]],
        title = metric,
        metric_name = metric,
        row_split = NULL,
        col_order = method_meta$order,
        legend_title = metric,
        method_meta = method_meta,
        cell_size_in = cell_size_in,
        show_box_axis = TRUE,
        show_legend = TRUE,
        show_significance = FALSE
      )
      if (is.null(ht)) {
        next
      }
      ht_dropout <- if (is.null(ht_dropout)) ht else (ht_dropout + ht)
    }
    if (!is.null(ht_dropout)) {
      heatmaps[["dropout"]] <- list(
        heatmap = ht_dropout,
        gap_in = gap_in
      )
    }


  }

  heatmaps
}

plot_heatmaps_combined <- function(
  synthetic_data,
  curated_data,
  show_significance = TRUE
) {
  heatmaps <- list()
  metrics <- intersect(
    resolve_plot_metrics(synthetic_data),
    resolve_plot_metrics(curated_data)
  )
  if (!length(metrics) || !isTRUE(has_complex_heatmap)) {
    return(invisible(NULL))
  }
  synthetic_data$Method <- normalize_beeline_method_name(synthetic_data$Method)
  curated_data$Method <- normalize_beeline_method_name(curated_data$Method)
  method_meta <- resolve_plot_methods(dplyr::bind_rows(
    synthetic_data,
    curated_data
  ))

  cell_med <- synthetic_cell_level_medians(synthetic_data, metrics)
  if (is.null(cell_med) || !nrow(cell_med)) {
    return(invisible(NULL))
  }
  cur <- curated_data
  cur$DatasetGroup <- resolve_curated_dataset_group(cur)

  families <- levels(cell_med$Family)
  groups <- sort(unique(as.character(cur$DatasetGroup)))

  cell_size_in <- 0.157
  gap_in <- 0.276

  mats_sub <- list()
  sub_splits <- list()
  row_order_sub <- c(families, groups)
  row_split_sub <- factor(
    c(rep("synthetic", length(families)), rep("curated", length(groups))),
    levels = c("synthetic", "curated")
  )

  for (metric in metrics) {
    syn_fam <- cell_med[cell_med$Metric == metric, , drop = FALSE] |>
      dplyr::group_by(Family, Method) |>
      dplyr::summarise(
        Value = stats::median(Value, na.rm = TRUE),
        .groups = "drop"
      )
    cur_grp <- cur |>
      dplyr::group_by(DatasetGroup, Method) |>
      dplyr::summarise(
        Value = stats::median(.data[[metric]], na.rm = TRUE),
        .groups = "drop"
      )
    syn_fam$RowKey <- as.character(syn_fam$Family)
    cur_grp$RowKey <- as.character(cur_grp$DatasetGroup)
    comb_df <- dplyr::bind_rows(
      syn_fam[, c("RowKey", "Method", "Value"), drop = FALSE],
      cur_grp[, c("RowKey", "Method", "Value"), drop = FALSE]
    )
    mat10 <- build_matrix_from_df(comb_df, "RowKey", "Method", "Value")
    mat10 <- mat10[row_order_sub, method_meta$order, drop = FALSE]
    mats_sub[[metric]] <- mat10
    sub_splits[[metric]] <- row_split_sub
  }

  ht_sub_list <- NULL
  for (i in seq_along(metrics)) {
    metric <- metrics[[i]]
    ht <- make_metric_heatmap(
      mat = mats_sub[[metric]],
      title = metric,
      metric_name = metric,
      row_split = sub_splits[[metric]],
      col_order = method_meta$order,
      legend_title = metric,
      method_meta = method_meta,
      cell_size_in = cell_size_in,
      show_box_axis = TRUE,
      show_legend = TRUE,
      show_significance = show_significance
    )
    if (is.null(ht)) {
      next
    }
    ht_sub_list <- if (is.null(ht_sub_list)) ht else (ht_sub_list + ht)
  }
  if (!is.null(ht_sub_list)) {
    heatmaps[["subdataset"]] <- list(
      heatmap = ht_sub_list,
      gap_in = gap_in
    )
  }

  heatmaps
}

resolve_plot_methods <- function(evaluation_data) {
  methods <- unique(as.character(evaluation_data$Method))
  method_order <- preferred_beeline_method_order(methods)
  list(
    methods = methods,
    order = method_order,
    colors = beeline_method_colors(method_order),
    shapes = beeline_method_shapes(method_order)
  )
}

resolve_plot_metrics <- function(evaluation_data) {
  metric_candidates <- c(
    "AUPRC",
    "AUROC",
    "ACC",
    "Precision",
    "Recall",
    "F1",
    "EarlyPrecision",
    "EPR"
  )
  excluded_metrics <- c("SI", "JI", "RuntimeSeconds", "PositiveRate")

  exclude <- c(
    "Dataset",
    "Method",
    "Case",
    "DatasetGroup",
    "ReferenceNetwork",
    "config_id",
    "cross_validation",
    "penalty",
    "r_threshold",
    "primary_metric_value",
    "SelectedMetric",
    excluded_metrics
  )
  available <- intersect(metric_candidates, colnames(evaluation_data))
  extra <- setdiff(colnames(evaluation_data), c(exclude, available))
  extra_numeric <- extra[vapply(
    extra,
    function(col) {
      vals <- suppressWarnings(as.numeric(evaluation_data[[col]]))
      mean(is.finite(vals)) > 0.9
    },
    logical(1)
  )]
  metrics <- c(available, extra_numeric)
  metrics <- metrics[metrics %in% colnames(evaluation_data)]
  unique(setdiff(metrics, excluded_metrics))
}

synthetic_cell_level_medians <- function(evaluation_data, metrics) {
  synthetic_info <- parse_synthetic_dataset_info(evaluation_data$Dataset)
  if (nrow(synthetic_info) != nrow(evaluation_data)) {
    return(NULL)
  }

  df <- cbind(evaluation_data, synthetic_info)
  df$CellsNum <- suppressWarnings(as.numeric(df$Cells))

  long <- tidyr::pivot_longer(
    df[, c("Method", "Family", "Cells", "CellsNum", metrics), drop = FALSE],
    cols = tidyselect::all_of(metrics),
    names_to = "Metric",
    values_to = "Value"
  )
  long$Value <- suppressWarnings(as.numeric(long$Value))
  long <- long[is.finite(long$Value), , drop = FALSE]

  med <- long |>
    dplyr::group_by(Method, Family, Cells, CellsNum, Metric) |>
    dplyr::summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

  fam_order <- sort(unique(as.character(med$Family)))
  cell_order <- sort(unique(med$CellsNum))
  med$Family <- factor(med$Family, levels = fam_order)
  med$Cells <- factor(
    med$Cells,
    levels = unique(as.character(med$Cells[order(med$CellsNum)]))
  )
  med
}

resolve_curated_dataset_group <- function(evaluation_data) {
  if ("DatasetGroup" %in% colnames(evaluation_data)) {
    return(evaluation_data$DatasetGroup)
  }
  parts <- strsplit(as.character(evaluation_data$Dataset), "-", fixed = TRUE)
  vapply(
    parts,
    function(p) if (length(p)) p[[1]] else NA_character_,
    character(1)
  )
}
