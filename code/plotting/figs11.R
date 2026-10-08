#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)
data_root <- file.path(
  repo,
  "results",
  "gse194122",
  "state_transitions",
  "transition_summary"
)
figure_root <- file.path(repo, "figures")
regulator_path <- file.path(
  data_root,
  "GSE5e_regulator_transition_counts_display.tsv"
)
assert_files(regulator_path)

regulators <- fread(regulator_path)
stopifnot(
  all(
    regulators$turnover ==
      regulators$gained + regulators$lost + regulators$shared_sign_flip
  ),
  all(regulators$net_edge_change == regulators$gained - regulators$lost)
)

transition_labels <- c(
  "HSC → MK/E",
  "MK/E → ProE",
  "ProE → Ery",
  "Ery → Norm"
)
regulators[, `:=`(
  transition_label = factor(
    transition_labels[transition_order],
    levels = transition_labels
  ),
  net_edge_change = as.numeric(net_edge_change),
  turnover = as.numeric(turnover)
)]

make_regulator_plot <- function(endpoint_name) {
  plot_data <- copy(regulators[endpoint == endpoint_name])
  regulator_order <- plot_data[,
    .(peak_turnover = max(turnover)),
    by = regulator
  ][order(peak_turnover, regulator), regulator]
  plot_data[, regulator := factor(regulator, levels = regulator_order)]
  max_abs_change <- max(abs(plot_data$net_edge_change))
  size_breaks <- pretty(c(0, max(plot_data$turnover)), n = 3)
  size_breaks <- size_breaks[
    size_breaks > 0 & size_breaks <= max(plot_data$turnover)
  ]
  if (length(size_breaks) > 3L) {
    size_breaks <- size_breaks[
      seq(1, length(size_breaks), length.out = 3)
    ]
  }

  ggplot(plot_data, aes(transition_label, regulator)) +
    geom_point(
      aes(size = turnover, fill = net_edge_change),
      shape = 21,
      stroke = 0.35,
      colour = INK_COLOUR
    ) +
    scale_size_continuous(
      name = "Turnover",
      range = c(0.8, 3.4),
      breaks = size_breaks
    ) +
    scale_fill_gradient2(
      name = "Gained − lost",
      low = "#3D75B5",
      mid = "white",
      high = "#DE615A",
      midpoint = 0,
      limits = c(-max_abs_change, max_abs_change),
      breaks = c(-max_abs_change, 0, max_abs_change),
      labels = label_number(big.mark = ",")
    ) +
    labs(x = NULL, y = NULL) +
    guides(
      size = guide_legend(order = 1, title.position = "top", nrow = 1),
      fill = guide_colourbar(
        display = "rectangles",
        order = 2,
        title.position = "top",
        barwidth = unit(0.45, "in"),
        barheight = unit(0.06, "in")
      )
    ) +
    theme_genome(6.8) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 5.9),
      axis.text.y = element_text(size = 5.6),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.spacing.x = unit(0.15, "in"),
      legend.title = element_text(size = 5.0),
      legend.text = element_text(size = 5.0),
      legend.key.size = unit(0.13, "in"),
      legend.margin = margin(t = 0, unit = "pt"),
      legend.box.margin = margin(0, 0, 0, -0.35, "in"),
      plot.margin = margin(3, 4, 1, 4)
    )
}

figs11a <- make_regulator_plot("TF-gene")
figs11b <- make_regulator_plot("factorized-chain")
save_pdf(figs11a, file.path(figure_root, "figs11a.pdf"), 1.5, 2.8)
save_pdf(figs11b, file.path(figure_root, "figs11b.pdf"), 1.5, 2.8)

source(file.path(repo, "code", "functions", "pdf_assembly.R"), local = TRUE)
assemble_pdf_figure(
  stem = "figs11",
  page = c(125, 115),
  panels = list(
    pdf_panel("a", 8, 111, 50),
    pdf_panel("b", 67, 111, 50)
  ),
  repo = repo,
  figure_root = figure_root
)
message("Saved Supplementary Fig. S11 complete-fit sensitivity panels.")
