#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(grid)
})
data.table::setDTthreads(1L)

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)
source(
  file.path(repo, "code", "functions", "beeline_plot_helpers.R"),
  local = TRUE
)

snapshot_root <- file.path(repo, "results", "beeline")
metric_root <- file.path(snapshot_root, "metrics")
runtime_root <- file.path(snapshot_root, "runtime")
source(file.path(repo, "code/functions/beeline_case_views.R"), local = TRUE)
figure_root <- file.path(repo, "figures")

paths <- c(
  metrics = file.path(
    metric_root,
    "case_metrics.tsv"
  ),
  ranks = file.path(metric_root, "rank_summary.tsv"),
  runtime = file.path(runtime_root, "runtime_all_methods.tsv"),
  runtime_summary = file.path(runtime_root, "runtime_summary.tsv")
)

comparison <- fread(paths[["metrics"]])
synthetic <- as.data.frame(comparison[DatasetType == "synthetic"])
curated <- as.data.frame(comparison[DatasetType == "curated"])
fig2a <- plot_heatmaps_combined(
  synthetic,
  curated,
  show_significance = FALSE
)$subdataset
if (is.null(fig2a)) {
  stop("Unable to construct Fig. 2a ComplexHeatmap.")
}
save_grid_pdf(
  function() {
    grid.newpage()
    ComplexHeatmap::draw(
      fig2a$heatmap,
      newpage = FALSE,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      ht_gap = unit(0.276, "in")
    )
  },
  file.path(figure_root, "fig2a.pdf"),
  8,
  5
)

cm <- copy(comparison)
cm[, rank_AUPRC := frank(-AUPRC, ties.method = "average"), by = Case]
cm[, rank_EPR := frank(-EPR, ties.method = "average"), by = Case]

ranks_long <- rbindlist(list(
  cm[, .(Case, Method, Metric = "AUPRC", Rank = rank_AUPRC)],
  cm[, .(Case, Method, Metric = "EPR", Rank = rank_EPR)]
))

ranks_long[, Best := Rank == min(Rank), by = .(Case, Metric)]

method_summary <- ranks_long[Metric == "AUPRC", .(
  median_rank = median(Rank),
  mean_rank = mean(Rank),
  top1_pct = mean(Best) * 100,
  top3_pct = mean(Rank <= 3 | Best) * 100
), by = Method][order(median_rank, mean_rank)]

method_order <- rev(method_summary$Method)
ranks_long[, Method := factor(Method, levels = method_order)]

stats <- ranks_long[, .(
  median = median(Rank),
  q25 = quantile(Rank, 0.25),
  q75 = quantile(Rank, 0.75),
  mean = mean(Rank),
  top1 = mean(Best) * 100,
  top3 = mean(Rank <= 3 | Best) * 100,
  cases = .N
), by = .(Metric, Method)]

colors <- method_colours[levels(ranks_long$Method)]

p_forest_auprc <- ggplot(stats[Metric == "AUPRC"], aes(y = Method)) +
  geom_segment(aes(x = q25, xend = q75, yend = Method, colour = Method), linewidth = 2.4, alpha = 0.85) +
  geom_point(aes(x = median, fill = Method), shape = 21, size = 3.0, colour = "black", stroke = 0.6) +
  geom_point(aes(x = mean), shape = 18, size = 2.0, colour = "#212121") +
  annotate(
    "text", x = 11, y = "multiCSN",
    label = "Median: 1.0",
    colour = "#D81B60", size = 2.1, fontface = "bold.italic", hjust = 1
  ) +
  scale_colour_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11), limits = c(0.8, 11.2)) +
  labs(title = "AUPRC rank", x = "Rank (1 = best)", y = NULL) +
  theme_genome(7.0) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#ECEFF1", linewidth = 0.3),
    panel.grid.major.y = element_line(colour = "#F5F5F5", linewidth = 0.3),
    plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 6.8),
    plot.margin = margin(2, 4, 2, 2)
  )

p_forest_epr <- ggplot(stats[Metric == "EPR"], aes(y = Method)) +
  geom_segment(aes(x = q25, xend = q75, yend = Method, colour = Method), linewidth = 2.4, alpha = 0.85) +
  geom_point(aes(x = median, fill = Method), shape = 21, size = 3.0, colour = "black", stroke = 0.6) +
  geom_point(aes(x = mean), shape = 18, size = 2.0, colour = "#212121") +
  annotate(
    "text", x = 11, y = "multiCSN",
    label = "Median: 1.5",
    colour = "#D81B60", size = 2.1, fontface = "bold.italic", hjust = 1
  ) +
  scale_colour_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11), limits = c(0.8, 11.2)) +
  labs(title = "EPR rank", x = "Rank (1 = best)", y = NULL) +
  theme_genome(7.0) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#ECEFF1", linewidth = 0.3),
    panel.grid.major.y = element_line(colour = "#F5F5F5", linewidth = 0.3),
    plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(2, 4, 2, 2)
  )

win_df <- stats[Metric == "AUPRC", .(
  Method,
  top1,
  top2_3 = top3 - top1,
  other = 100 - top3
)]
win_melt <- data.table::melt(win_df, id.vars = "Method", variable.name = "Tier", value.name = "Percent")
win_melt[, Tier := factor(Tier, levels = c("other", "top2_3", "top1"))]

p_win <- ggplot(win_melt, aes(y = Method, x = Percent, fill = Tier)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.25) +
  annotate(
    "text", x = 26, y = "multiCSN",
    label = sprintf("%.1f%% Best", win_df[as.character(Method) == "multiCSN", top1]), colour = "white", size = 2.0, fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("top1" = "#D81B60", "top2_3" = "#42A5F5", "other" = "#E0E0E0"),
    labels = c("top1" = "Best (incl. ties)", "top2_3" = "Rank 2\u20133", "other" = "Rank \u22654"),
    name = NULL
  ) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100), breaks = c(0, 50, 100)) +
  labs(title = "Top-tier frequency", x = "% of cases", y = NULL) +
  theme_genome(7.0) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 5.8),
    legend.margin = margin(l = -4, r = 2),
    legend.key.size = unit(0.08, "in"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#ECEFF1", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(2, 2, 2, 2)
  )

fig2b <- (p_forest_auprc | p_forest_epr | p_win) +
  plot_layout(widths = c(1.25, 1.05, 1.25))
save_pdf(fig2b, file.path(figure_root, "fig2b.pdf"), 8, 2.0)

source(file.path(repo, "code/plotting/fig2c.R"), local = TRUE)

case_view <- read_beeline_case(repo, "HSC-2000-3-30")
truth <- case_view$truth
edges <- case_view$edges
summary <- case_view$summary
display_methods <- c(
  "multiCSN",
  "GENIE3",
  "GRNBoost2",
  "PIDC",
  "PPCOR",
  "SCODE",
  "SCORPION",
  "SINCERITIES",
  "LEAP",
  "GNIPLR",
  "SCSGL"
)
budget <- edges[Method == "multiCSN", .N]
regulators <- unique(truth$regulator)
targets <- unique(truth$target)
truth_table <- truth[, .(regulator, target, weight = as.numeric(sign))]
selected_edges <- setNames(lapply(display_methods, function(method) {
  edges[Method == method]
}), display_methods)
network_tables <- lapply(selected_edges, function(selected) {
  selected[, .(
    regulator,
    target,
    weight = as.numeric(within_network_signed_rank)
  )]
})
summary_lookup <- setNames(summary$AUPRC, summary$Method)
heatmap_titles <- c(
  "Ground truth",
  sprintf("%s (AUPRC %.3f)", display_methods, summary_lookup[display_methods])
)
fig2d <- multiCSN::plot_network_heatmap(
  c(list("Ground truth" = truth_table), network_tables),
  regulators = regulators,
  targets = targets,
  show_names = TRUE,
  show_names_position = "all",
  rect_color = "grey92",
  heatmap_size_lock = TRUE,
  performance_metrics = character(),
  truth_cell_border_lwd = 2,
  nrow = 3,
  heatmap_title = heatmap_titles,
  heatmap_title_color = c("black", unname(method_colours[display_methods])),
  heatmap_palcolor = c("#2166AC", "#F7F7F7", "#B2182B"),
  legend_name = c("Truth", rep("Weight", length(network_tables))),
  row_title = "Regulators"
)
save_grid_pdf(
  function() {
    grid.newpage()
    grid.draw(italicise_gene_grobs(grid::forceGrob(fig2d), c(regulators, targets)))
  },
  file.path(figure_root, "fig2d.pdf"),
  15,
  9
)

source(file.path(repo, "code", "functions", "pdf_assembly.R"), local = TRUE)
assemble_pdf_figure(
  stem = "fig2",
  page = c(183, 225),
  panels = list(
    pdf_panel("a", 8, 221, 111),
    pdf_panel("b", 8, 160.5, 171),
    pdf_panel("c", 121, 221, 58),
    pdf_panel("d", 8, 118.0, 171)
  ),
  repo = repo,
  figure_root = figure_root
)
