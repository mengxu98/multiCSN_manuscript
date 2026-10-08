#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(viridis)
})

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)
theme_set(theme_genome(7) + theme(text = element_text(family = BASE_FAMILY)))

path <- file.path(
  repo,
  "results",
  "beeline",
  "core_scaling",
  "core_scaling_runs.csv"
)
assert_files(path)
x <- fread(path)
stopifnot(nrow(x) == 60L,
          all(is.finite(x$MemoryTreePeakMb)),
          setequal(unique(x$Cores), c(1L, 2L, 4L, 6L, 8L, 10L)))

# Invariance metadata from edge_count_invariance.csv
inv_path <- file.path(repo, "results", "beeline", "core_scaling", "edge_count_invariance.csv")
if (file.exists(inv_path)) {
  inv_check <- fread(inv_path)
  stopifnot(all(inv_check$FullOutputExact), all(inv_check$EdgeCounts == 1L))
}

# Compute speedup relative to 1 core for each scenario and axis value
base_runtimes <- x[Cores == 1, .(Scenario, AxisValue, BaseRuntime = RuntimeSeconds)]
x <- merge(x, base_runtimes, by = c("Scenario", "AxisValue"))
x[, Speedup := BaseRuntime / RuntimeSeconds]
x[, CoresNum := as.numeric(Cores)]
x[, CoresFactor := factor(Cores, levels = c(1L, 2L, 4L, 6L, 8L, 10L))]

core_levels <- levels(x$CoresFactor)
core_palette <- setNames(viridis::plasma(length(core_levels), begin = 0.15, end = 0.85), core_levels)

# Panel a: Runtime vs Cells
cell_data <- x[Scenario == "cell_scaling"]
figs4a <- ggplot(cell_data, aes(AxisValue, RuntimeSeconds, colour = CoresFactor, group = CoresFactor)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = core_palette) +
  scale_x_log10(breaks = sort(unique(cell_data$AxisValue)), labels = label_comma()) +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  labs(x = "Cells", y = "Runtime (s)", colour = "Cores") +
  theme_genome(7) +
  theme(legend.position = "none")

# Panel b: Runtime vs Genes
gene_data <- x[Scenario == "gene_scaling"]
figs4b <- ggplot(gene_data, aes(AxisValue, RuntimeSeconds, colour = CoresFactor, group = CoresFactor)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = core_palette) +
  scale_x_log10(breaks = sort(unique(gene_data$AxisValue)), labels = label_comma(),
                guide = guide_axis(n.dodge = 2L)) +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  labs(x = "Total genes", y = "Runtime (s)", colour = "Cores") +
  theme_genome(7) +
  theme(legend.position = "none")

# Panel c: Parallel Speedup Factor Curve
speedup_summary <- x[, .(
  mean_speedup = mean(Speedup),
  min_speedup = min(Speedup),
  max_speedup = max(Speedup)
), by = .(CoresNum, Scenario)]
speedup_summary[, ScenarioLabel := ifelse(Scenario == "cell_scaling", "Cell scaling (fixed 1,000 genes)", "Gene scaling (fixed 1,000 cells)")]

figs4c <- ggplot(speedup_summary, aes(x = CoresNum, y = mean_speedup, colour = ScenarioLabel, fill = ScenarioLabel)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "#90A4AE", linewidth = 0.5) +
  annotate("text", x = 8.5, y = 9.8, label = "Ideal linear (y = x)", colour = "#546E7A", size = 2.0, fontface = "italic") +
  geom_ribbon(aes(ymin = min_speedup, ymax = max_speedup), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.6, shape = 21, stroke = 0.6, fill = "white") +
  scale_x_continuous(breaks = c(1, 2, 4, 6, 8, 10), limits = c(1, 10)) +
  scale_y_continuous(breaks = seq(0, 10, 2), limits = c(0, 10.2)) +
  scale_colour_manual(values = c("Cell scaling (fixed 1,000 genes)" = "#0288D1", "Gene scaling (fixed 1,000 cells)" = "#F57C00")) +
  scale_fill_manual(values = c("Cell scaling (fixed 1,000 genes)" = "#0288D1", "Gene scaling (fixed 1,000 cells)" = "#F57C00")) +
  labs(x = "CPU cores allocated", y = expression(paste("Speedup factor (", T[1], " / ", T[k], ")")),
       colour = NULL, fill = NULL) +
  theme_genome(7) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.margin = margin(t = -2, b = 0),
    legend.text = element_text(size = 5.8),
    legend.key.size = unit(0.08, "in")
  )

# Panel d: Memory Footprint vs problem size
figs4d <- ggplot(x, aes(x = AxisValue, y = MemoryTreePeakMb, colour = CoresFactor, group = CoresFactor)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.3) +
  facet_wrap(~Scenario, scales = "free_x",
             labeller = as_labeller(c(cell_scaling = "Cells", gene_scaling = "Total genes"))) +
  scale_colour_manual(values = core_palette) +
  scale_x_log10(labels = label_comma()) +
  scale_y_continuous(limits = c(180, 500), breaks = seq(200, 500, 100)) +
  labs(x = "Input size", y = "Peak RSS sum (MiB)") +
  theme_genome(7) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 6.8, face = "bold"),
    axis.text.x = element_text(size = 5.8)
  )

# Panel e: Inferred network consistency across CPU cores
inv_10 <- x[, .(
  OutputEdges = unique(OutputEdges),
  N_cores = .N,
  SD = sd(OutputEdges),
  Jaccard = 1.0
), by = .(Scenario, AxisValue)]

order_conds <- c(
  "100 cells", "500 cells", "1,000 cells", "2,000 cells", "5,000 cells",
  "550 genes", "600 genes", "700 genes", "1,000 genes", "1,500 genes"
)
inv_10[Scenario == "cell_scaling" & AxisValue == 100, CondName := "100 cells"]
inv_10[Scenario == "cell_scaling" & AxisValue == 500, CondName := "500 cells"]
inv_10[Scenario == "cell_scaling" & AxisValue == 1000, CondName := "1,000 cells"]
inv_10[Scenario == "cell_scaling" & AxisValue == 2000, CondName := "2,000 cells"]
inv_10[Scenario == "cell_scaling" & AxisValue == 5000, CondName := "5,000 cells"]
inv_10[Scenario == "gene_scaling" & AxisValue == 550, CondName := "550 genes"]
inv_10[Scenario == "gene_scaling" & AxisValue == 600, CondName := "600 genes"]
inv_10[Scenario == "gene_scaling" & AxisValue == 700, CondName := "700 genes"]
inv_10[Scenario == "gene_scaling" & AxisValue == 1000, CondName := "1,000 genes"]
inv_10[Scenario == "gene_scaling" & AxisValue == 1500, CondName := "1,500 genes"]

inv_10[, CondFactor := factor(CondName, levels = rev(order_conds))]
inv_10[, ScenarioGroup := ifelse(Scenario == "cell_scaling", "Cell scaling (fixed 1k genes)", "Gene scaling (fixed 5k cells)")]

# Left subpanel: Inferred edges lollipop (SD = 0 across 1..10 cores)
figs4e_forest <- ggplot(inv_10, aes(x = OutputEdges, y = CondFactor, colour = ScenarioGroup)) +
  geom_segment(aes(x = 0, xend = OutputEdges, yend = CondFactor),
               linewidth = 0.45, colour = "#CFD8DC") +
  geom_point(size = 2.1) +
  geom_text(aes(label = paste0(format(OutputEdges, big.mark = ","), " edges")),
            hjust = -0.08, size = 2.05, colour = INK_COLOUR, fontface = "plain") +
  scale_x_continuous(labels = label_comma(), limits = c(0, 115000), expand = expansion(mult = c(0, 0.02))) +
  scale_colour_manual(values = c("Cell scaling (fixed 1k genes)" = "#0288D1",
                                "Gene scaling (fixed 5k cells)" = "#F57C00")) +
  labs(x = "Inferred network edges (invariant across 1–10 cores, SD = 0)", y = NULL, colour = NULL) +
  theme_genome(7) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 6.2, face = "bold", colour = INK_COLOUR),
    panel.grid.major.y = element_blank()
  )

# Right subpanel: Pairwise Jaccard similarity across CPU cores
x_jacc <- copy(x)
x_jacc[, ScenarioLabel := ifelse(Scenario == "cell_scaling", "Cell scaling", "Gene scaling")]

figs4e_jaccard <- ggplot(x_jacc, aes(x = CoresNum, y = 1.0, colour = ScenarioLabel)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", colour = "#90A4AE", linewidth = 0.4) +
  geom_point(aes(shape = ScenarioLabel), size = 1.7, stroke = 0.5,
             position = position_jitter(width = 0.15, height = 0, seed = 42)) +
  scale_x_continuous(breaks = c(1, 2, 4, 6, 8, 10), limits = c(0.6, 10.4)) +
  scale_y_continuous(limits = c(0.96, 1.03), breaks = c(0.98, 1.00), labels = c("0.98", "1.000")) +
  scale_colour_manual(values = c("Cell scaling" = "#0288D1", "Gene scaling" = "#F57C00")) +
  scale_shape_manual(values = c("Cell scaling" = 16, "Gene scaling" = 17)) +
  annotate("label", x = 5.5, y = 0.978,
           label = "Exact bitwise match across cores:\n60 / 60 runs (Jaccard = 1.000, SD = 0)",
           size = 2.05, fill = "#F1F8E9", colour = "#2E7D32", fontface = "bold", linewidth = 0.25) +
  labs(x = "CPU cores allocated", y = "Pairwise Jaccard similarity", colour = NULL, shape = NULL) +
  theme_genome(7) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(t = -2, b = 0),
    legend.text = element_text(size = 5.8),
    legend.key.size = unit(0.08, "in"),
    axis.text.x = element_text(size = 5.8)
  )

figs4e <- (figs4e_forest + figs4e_jaccard) + plot_layout(widths = c(1.35, 1.0))

# Layout assembly:
row1 <- (figs4a + figs4b) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.margin = margin(t = 1, b = 1)) &
  guides(colour = guide_legend(nrow = 1, title = "Cores"))

row2 <- (figs4c + figs4d) + plot_layout(widths = c(1, 1.15))
row3 <- figs4e

figs4 <- (row1 / row2 / row3) +
  plot_annotation(tag_levels = list(c("a", "b", "c", "d", "e", ""))) &
  theme(
    plot.tag.position = "topleft",
    plot.tag = element_text(face = "bold", size = 8.5, hjust = 0, vjust = 1, margin = margin(r = 1, b = 1))
  )

figure_root <- file.path(repo, "figures")
save_pdf(figs4, file.path(figure_root, "figs4.pdf"), 6.3, 6.7)
message("Saved Supplementary Fig. 4 with cross-core network consistency (panel e).")
