#!/usr/bin/env Rscript
# Refresh only the runtime panel; reuse the existing accuracy and topology PDFs.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
repo <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(repo, "code/functions/common.R"), local = TRUE)
runtime_root <- file.path(repo, "results/beeline/runtime")
figure_root <- file.path(repo, "figures")
assert_files(file.path(runtime_root, "runtime_provenance.tsv"))
runtime <- fread(file.path(runtime_root, "runtime_all_methods.tsv"))
runtime_summary <- fread(file.path(runtime_root, "runtime_summary.tsv"))
method_order <- c(
  "multiCSN",
  runtime_summary[Method != "multiCSN"][order(total_runtime_seconds)]$Method
)
runtime[, Method := factor(Method, levels = rev(method_order))]
runtime_summary[, Method := factor(Method, levels = rev(method_order))]
fig2c_left <- ggplot(
  runtime_summary,
  aes(total_runtime_seconds, Method, colour = Method)
) +
  geom_segment(aes(x = 1, xend = total_runtime_seconds, yend = Method),
               linewidth = 1.6, alpha = 0.88) +
  geom_point(size = 1.5) +
  geom_text(
    aes(x = total_runtime_seconds * 1.18,
        label = sprintf("%.1fx", runtime_ratio_to_multicsn)),
    hjust = 0,
    size = 6.0 / ggplot2::.pt,
    colour = "grey25"
  ) +
  scale_colour_manual(values = method_colours, drop = FALSE) +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000, 10000),
    labels = c("1", "10", "100", "1,000", "10,000"),
    expand = expansion(mult = c(0.02, 0.30))
  ) +
  labs(x = "Total runtime (s)", y = NULL) +
  theme_genome(6.8) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 6.2),
    axis.text.x = element_text(size = 6.0),
    axis.title.x = element_text(size = 6.8),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(1.5, 2, 1.5, 1.5)
  )
fig2c_right <- ggplot(runtime, aes(RuntimeSeconds, Method, fill = Method)) +
  geom_boxplot(
    width = 0.65,
    alpha = 0.88,
    outlier.size = 0.20,
    outlier.alpha = 0.25,
    linewidth = 0.20
  ) +
  scale_fill_manual(values = method_colours, drop = FALSE) +
  scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
  labs(x = "Runtime per case (s)", y = NULL) +
  theme_genome(6.8) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 6.2),
    axis.text.x = element_text(size = 6.0),
    axis.title.x = element_text(size = 6.8),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(1.5, 2, 1.5, 1.5)
  )
fig2c <- fig2c_left / fig2c_right + plot_layout(heights = c(1, 1))
save_pdf(fig2c, file.path(figure_root, "fig2c.pdf"), 2.55, 2.60)

if (sys.nframe() == 0L) {
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
}
