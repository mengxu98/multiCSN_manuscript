suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
  library(Matrix)
  library(SeuratObject)
  library(scop)
})

class_order <- c("GluN", "IN", "IPC_RG", "OPC_Glia", "Vascular")
class_labels <- c("GluN", "IN", "IPC/RG", "OPC/Glia", "Vascular")
class_colours <- c(
  GluN = "#3B6FB6", IN = "#8A5BA7", IPC_RG = "#D98A2B",
  OPC_Glia = "#4AA684", Vascular = "#C4474D"
)

region_label <- function(region) {
  x <- tstrsplit(sub("^chr", "", region), "-", fixed = TRUE)
  paste0(
    "chr", x[[1L]], ":", comma(as.numeric(x[[2L]])), "–",
    comma(as.numeric(x[[3L]]))
  )
}

read_locus_inputs <- function(repo, target) {
  base <- file.path(repo, "results", "gse162170", "loci", "candidates")
  files <- setNames(
    file.path(base, paste0(target, c(
      "_accessibility.tsv", "_expression_summary.tsv",
      "_expression_samples.tsv", "_expression_cells.tsv",
      "_gene_annotation.tsv", "_paths.tsv"
    ))),
    c(
      "accessibility", "expression", "expression_samples",
      "expression_cells", "annotation", "paths"
    )
  )
  assert_files(files)
  c(list(target = target), lapply(files, fread))
}

locus_limits <- function(obj) {
  range(c(obj$paths$start, obj$paths$end, obj$annotation$start, obj$annotation$end))
}

make_arc_data <- function(paths, target_tss) {
  paths <- copy(paths)
  paths[, midpoint := (start + end) / 2]
  setorder(paths, start, end, -chain_delta_bic, regulator)
  paths[, path_id := seq_len(.N)]
  paths[, arc_height := 0.35 + 0.08 * frank(
    chain_delta_bic,
    ties.method = "first", na.last = "keep"
  )]
  paths[, sign_label := ifelse(chain_direction >= 0, "Positive", "Negative")]
  rbindlist(lapply(seq_len(nrow(paths)), function(i) {
    t <- seq(0, 1, length.out = 100L)
    data.table(
      path_id = paths$path_id[[i]],
      x = (1 - t) * paths$midpoint[[i]] + t * target_tss,
      y = 4 * paths$arc_height[[i]] * t * (1 - t),
      sign_label = paths$sign_label[[i]],
      chain_delta_bic = paths$chain_delta_bic[[i]]
    )
  }))
}

make_locus_paths <- function(obj) {
  paths <- copy(obj$paths)
  annotation <- copy(obj$annotation)
  strand <- unique(annotation$strand[annotation$strand %chin% c("+", "-")])
  stopifnot(length(strand) == 1L)
  target_tss <- if (strand == "+") min(annotation$start) else max(annotation$end)
  curves <- make_arc_data(paths, target_tss)
  sign_colours <- c(Negative = "#3B86B6", Positive = MULTICSN_COLOUR)

  ggplot(curves, aes(x, y,
    group = path_id, colour = sign_label,
    linewidth = chain_delta_bic
  )) +
    geom_path(alpha = 0.72, lineend = "round") +
    geom_point(
      data = unique(paths[, .(midpoint = (start + end) / 2)]),
      aes(midpoint, 0), inherit.aes = FALSE, shape = 23, size = 1.55,
      fill = "#62B44A", colour = "#2E6C29", stroke = 0.35
    ) +
    geom_point(
      data = data.table(target_tss), aes(target_tss, 0), inherit.aes = FALSE,
      shape = 21, size = 1.8, fill = "#E85B24", colour = "#99310F", stroke = 0.35
    ) +
    scale_colour_manual(values = sign_colours, name = "Direction") +
    scale_linewidth_continuous(
      range = c(0.25, 0.9), trans = "log1p",
      name = expression(Delta * BIC)
    ) +
    guides(
      colour = guide_legend(
        ncol = 1, byrow = TRUE, order = 1,
        title.position = "top", title.hjust = 0
      ),
      linewidth = guide_legend(
        ncol = 1, byrow = TRUE, order = 2,
        title.position = "top", title.hjust = 0
      )
    ) +
    scale_x_continuous(
      labels = label_number(scale = 1e-6, accuracy = 0.01),
      limits = locus_limits(obj), expand = c(0, 0)
    ) +
    scale_y_continuous(NULL, breaks = NULL, expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = NULL) +
    theme_genome(7) +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
      panel.grid = element_blank(), legend.position = "right",
      legend.direction = "vertical", legend.box = "vertical",
      legend.box.just = "left", legend.key.width = unit(0.22, "in"),
      legend.spacing.y = unit(0.04, "in"), plot.title = element_blank(),
      plot.subtitle = element_blank(), plot.margin = margin(2, 4, 0, 4)
    )
}

make_locus_accessibility <- function(obj) {
  x <- copy(obj$accessibility)
  x[, broad_class := factor(broad_class, levels = class_order)]
  x[, class_index := match(as.character(broad_class), class_order)]
  maximum <- max(x$mean_normalized_signal, na.rm = TRUE)
  x[, scaled := if (maximum > 0) mean_normalized_signal / maximum else 0]
  baseline <- data.table(class_index = seq_along(class_order))
  limits <- locus_limits(obj)

  ggplot() +
    geom_segment(
      data = baseline,
      aes(x = limits[[1L]], xend = limits[[2L]], y = class_index, yend = class_index),
      linewidth = 0.22, colour = "#C8CDD3"
    ) +
    geom_segment(
      data = x,
      aes(
        x = (start + end) / 2, xend = (start + end) / 2,
        y = class_index, yend = class_index + 0.72 * scaled,
        colour = broad_class
      ),
      linewidth = 0.95, alpha = 0.88, lineend = "round"
    ) +
    scale_colour_manual(values = class_colours, guide = "none") +
    scale_y_continuous(
      breaks = seq_along(class_order), labels = class_labels,
      limits = c(0.72, length(class_order) + 0.92), expand = c(0, 0)
    ) +
    scale_x_continuous(
      labels = label_number(scale = 1e-6, accuracy = 0.01),
      limits = limits, expand = c(0, 0)
    ) +
    labs(x = NULL, y = NULL) +
    theme_genome(7) +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
      panel.grid = element_blank(), plot.title = element_blank(),
      plot.subtitle = element_blank(), plot.margin = margin(0, 4, 0, 4)
    )
}

make_locus_gene <- function(obj) {
  annotation <- copy(obj$annotation)
  exons <- unique(annotation[type == "exon", .(start, end)])
  gene_start <- min(annotation$start)
  gene_end <- max(annotation$end)
  strand <- unique(annotation$strand[annotation$strand %chin% c("+", "-")])
  regions <- unique(obj$paths[, .(start, end)])
  arrow_spec <- if (strand == "+") {
    arrow(type = "closed", length = unit(0.047, "in"))
  } else {
    arrow(type = "closed", length = unit(0.047, "in"), ends = "first")
  }

  ggplot() +
    geom_rect(
      data = regions, aes(xmin = start, xmax = end, ymin = 0.72, ymax = 0.94),
      fill = "#B6BEC8", colour = "#89929D", linewidth = 0.15
    ) +
    geom_segment(
      data = regions, aes(
        x = (start + end) / 2, xend = (start + end) / 2,
        y = 0.69, yend = 0.97
      ),
      colour = "#89929D", linewidth = 0.48
    ) +
    geom_segment(aes(x = gene_start, xend = gene_end, y = 0.34, yend = 0.34),
      linewidth = 0.5, colour = "#273A70", arrow = arrow_spec
    ) +
    geom_rect(
      data = exons, aes(xmin = start, xmax = end, ymin = 0.21, ymax = 0.47),
      fill = "#3C62A6", colour = "#273A70", linewidth = 0.12
    ) +
    annotate("text",
      x = (gene_start + gene_end) / 2, y = 0.12,
      label = paste0(obj$target, " (", strand, ")"), family = BASE_FAMILY,
      size = 6.5 / ggplot2::.pt, colour = "#273A70"
    ) +
    scale_x_continuous(
      labels = label_number(scale = 1e-6, accuracy = 0.01),
      limits = locus_limits(obj), expand = c(0, 0)
    ) +
    scale_y_continuous(NULL, breaks = NULL, limits = c(0, 1.08), expand = c(0, 0)) +
    labs(x = paste0(unique(obj$paths$chromosome), " position (Mb)")) +
    theme_genome(7) +
    theme(
      axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      panel.grid = element_blank(), plot.title = element_blank(),
      plot.subtitle = element_blank(), plot.margin = margin(0, 4, 2, 4)
    )
}

make_locus_track <- function(obj) {
  make_locus_paths(obj) / make_locus_accessibility(obj) / make_locus_gene(obj) +
    plot_layout(heights = c(1.6, 0.95, 0.88), guides = "collect")
}

make_locus_key <- function(obj) {
  x <- copy(obj$paths)
  setorder(x, start, end, -chain_delta_bic, regulator)
  stopifnot(all(is.finite(x$chain_delta_bic)), all(x$chain_delta_bic >= 0))
  x[, sign_label := ifelse(chain_direction >= 0, "Positive", "Negative")]
  x[, sign := ifelse(sign_label == "Positive", "+", "−")]
  x[, delta_label := fifelse(
    chain_delta_bic > 0 & chain_delta_bic < 0.01,
    "<0.01", sprintf("%.2f", chain_delta_bic)
  )]
  blocks <- if (nrow(x) > 20L) 2L else 1L
  region_groups <- x[, .(
    group_rows = .N,
    region_start = min(start)
  ), by = region][order(region_start)]
  if (blocks == 2L) {
    split_at <- which.min(abs(cumsum(region_groups$group_rows) - nrow(x) / 2))
    region_groups[, block := fifelse(seq_len(.N) <= split_at, 1L, 2L)]
  } else {
    region_groups[, block := 1L]
  }
  x[region_groups, block := i.block, on = "region"]

  block_layout <- lapply(seq_len(blocks), function(block_id) {
    block_regions <- region_groups[block == block_id]
    row_list <- vector("list", nrow(block_regions))
    region_list <- vector("list", nrow(block_regions))
    cursor <- 0
    for (i in seq_len(nrow(block_regions))) {
      region_id <- block_regions$region[[i]]
      region_rows <- copy(x[region == region_id])
      region_rows[, y := -(cursor + seq_len(.N) - 0.5)]
      row_list[[i]] <- region_rows
      region_list[[i]] <- data.table(
        block = block_id,
        region = region_id,
        y = mean(region_rows$y)
      )
      cursor <- cursor + nrow(region_rows) + 0.80
    }
    list(
      rows = rbindlist(row_list),
      regions = rbindlist(region_list),
      span = cursor - 0.80
    )
  })
  x <- rbindlist(lapply(block_layout, `[[`, "rows"))
  region_centres <- rbindlist(lapply(block_layout, `[[`, "regions"))
  max_span <- max(vapply(block_layout, `[[`, numeric(1), "span"))

  block_step <- 0.455
  x[, x_base := (block - 1L) * block_step]
  x[, `:=`(
    tf_x = x_base + 0.010,
    line_x = x_base + pmin(0.100, 0.017 + 0.0065 * nchar(regulator)),
    chr_node_x = x_base + 0.142,
    bic_x = x_base + 0.385,
    sign_x = x_base + 0.420
  )]
  region_centres[, x_base := (block - 1L) * block_step]
  region_centres[, `:=`(
    chr_node_x = x_base + 0.142,
    chr_text_x = x_base + 0.153,
    region_label = vapply(region, region_label, character(1))
  )]
  x[region_centres, region_y := i.y, on = c("block", "region")]
  column_headers <- rbindlist(lapply(seq_len(blocks), function(block_id) {
    x_base <- (block_id - 1L) * block_step
    data.table(
      x = c(x_base + 0.010, x_base + 0.153, x_base + 0.385, x_base + 0.420),
      y = 0.82,
      label = c("TF", "Chr", "ΔBIC", "Sign"),
      hjust = c(0, 0, 1, 0.5)
    )
  }))
  sign_colours <- c(Negative = "#3B86B6", Positive = MULTICSN_COLOUR)

  ggplot(x) +
    geom_curve(
      aes(x = line_x, xend = chr_node_x, y = y, yend = region_y),
      colour = "#B9C2C8",
      linewidth = 0.35,
      curvature = 0.16,
      angle = 90,
      lineend = "round"
    ) +
    geom_point(
      data = region_centres,
      aes(chr_node_x, y),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.9,
      stroke = 0.35,
      fill = "white",
      colour = "#7B878F"
    ) +
    geom_text(
      data = region_centres,
      aes(chr_text_x, y, label = region_label),
      inherit.aes = FALSE,
      hjust = 0,
      family = BASE_FAMILY,
      colour = INK_COLOUR,
      size = 5.2 / ggplot2::.pt
    ) +
    geom_text(
      aes(tf_x, y, label = regulator),
      hjust = 0,
      family = BASE_FAMILY,
      colour = INK_COLOUR,
      size = 5.7 / ggplot2::.pt
    ) +
    geom_text(
      aes(bic_x, y, label = delta_label),
      hjust = 1,
      family = BASE_FAMILY,
      colour = INK_COLOUR,
      size = 5.5 / ggplot2::.pt
    ) +
    geom_point(
      aes(sign_x, y, fill = sign_label),
      shape = 21,
      size = 2.35,
      stroke = 0.25,
      colour = "white"
    ) +
    geom_text(
      aes(sign_x, y, label = sign),
      family = BASE_FAMILY,
      fontface = "plain",
      colour = "white",
      size = 5.0 / ggplot2::.pt
    ) +
    geom_text(
      data = column_headers,
      aes(x, y, label = label, hjust = hjust),
      inherit.aes = FALSE,
      family = BASE_FAMILY,
      fontface = "plain",
      colour = INK_COLOUR,
      size = 5.8 / ggplot2::.pt
    ) +
    geom_segment(
      data = if (blocks == 2L) {
        data.table(
          x = block_step - 0.005, xend = block_step - 0.005,
          y = 1.08, yend = -max_span - 0.2
        )
      } else {
        data.table()
      },
      aes(x, y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      colour = "#DDE2E5",
      linewidth = 0.3
    ) +
    scale_fill_manual(values = sign_colours, guide = "none") +
    scale_x_continuous(limits = c(0, 0.895), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-max_span - 0.35, 1.20), expand = c(0, 0)) +
    theme_void(base_family = BASE_FAMILY, base_size = 7) +
    theme(
      plot.margin = margin(2, 3, 2, 3), plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
}

make_locus_expression <- function(obj) {
  x <- copy(obj$expression_samples)
  x[, cell_type := factor(broad_class, levels = class_order, labels = class_labels)]
  stopifnot(nrow(x) == 15L, uniqueN(x$sample_id) == 3L)
  ggplot(x, aes(cell_type, mean_normalized_signal, group = sample_id)) +
    geom_line(colour = "#B8C0C6", linewidth = 0.35) +
    geom_point(aes(fill = cell_type), shape = 21, size = 2.2,
      colour = "white", stroke = 0.3) +
    scale_fill_manual(values = setNames(class_colours, class_labels), guide = "none") +
    labs(x = NULL, y = "Sample-mean RNA expression") +
    theme_genome(7) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid = element_blank(), plot.title = element_blank(),
      plot.subtitle = element_blank(), legend.position = "none"
    )
}
