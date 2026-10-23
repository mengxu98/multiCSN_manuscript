#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)

# Explicitly set session theme for self-contained execution.
theme_set(theme_genome(7) + theme(text = element_text(family = BASE_FAMILY)))
endpoint_path <- file.path(
  repo,
  "results", "scmultisim", "multiseed", "per_seed_diagnostics.tsv"
)
component_path <- file.path(
  repo,
  "results", "scmultisim", "component_ablation",
  "endpoint_metrics_all_seeds.tsv"
)

x <- fread(endpoint_path)

scope_contract <- c(
  tf_static = "TF-gene:static_input_GRN",
  tf_dynamic = "TF-gene:dynamic_union",
  region_gene = "region-gene:all_variable_regions"
)

endpoint_colours <- c(
  `TF–gene` = "#2F6FA5",
  `Region–gene` = "#C47A2C",
  `Dynamic-union truth` = "#4B9C8E"
)

static_region <- x[
  scope %chin% scope_contract[c("tf_static", "region_gene")]
][,
  endpoint := fifelse(
    scope == scope_contract[["tf_static"]],
    "TF–gene",
    "Region–gene"
  )
]
static_region[,
  endpoint := factor(endpoint, levels = c("TF–gene", "Region–gene"))
]

endpoint_summary <- static_region[,
  .(
    prevalence = median(prevalence),
    average_precision = median(average_precision),
    average_precision_q25 = quantile(average_precision, 0.25, names = FALSE),
    average_precision_q75 = quantile(average_precision, 0.75, names = FALSE),
    auroc = median(auroc),
    auroc_q25 = quantile(auroc, 0.25, names = FALSE),
    auroc_q75 = quantile(auroc, 0.75, names = FALSE),
    ap_over_prevalence = median(average_precision / prevalence)
  ),
  by = endpoint
]

label_positions <- data.table(
  endpoint = factor(
    c("TF–gene", "Region–gene"),
    levels = levels(static_region$endpoint)
  ),
  label_x = c(0.0055, 0.0035),
  label_y = c(0.85, 0.12)
)
endpoint_summary <- merge(endpoint_summary, label_positions, by = "endpoint")
endpoint_summary[,
  label := sprintf(
    "%s\nmedian AP %.3f\n%.2f× baseline",
    endpoint,
    average_precision,
    ap_over_prevalence
  )
]

panel_title_theme <- theme(
  text = element_text(family = "Arial"),
  plot.title = element_text(
    family = "Arial",
    size = 7.2,
    face = "plain",
    colour = INK_COLOUR,
    margin = margin(0, 0, 4, 0)
  ),
  plot.tag = element_text(size = 8.5, face = "bold", colour = "black")
)

figs6a <- ggplot(
  static_region,
  aes(prevalence, average_precision, colour = endpoint)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.4
  ) +
  geom_segment(
    aes(xend = prevalence, y = prevalence, yend = average_precision),
    linewidth = 0.35,
    alpha = 0.25,
    show.legend = FALSE
  ) +
  geom_point(
    aes(y = prevalence),
    shape = 21,
    size = 1.5,
    stroke = 0.4,
    fill = "white",
    alpha = 0.8,
    show.legend = FALSE
  ) +
  geom_point(size = 1.65, alpha = 0.82, show.legend = FALSE) +
  geom_point(
    data = endpoint_summary,
    aes(prevalence, average_precision, fill = endpoint),
    inherit.aes = FALSE,
    shape = 21,
    size = 3.1,
    stroke = 0.75,
    colour = INK_COLOUR,
    show.legend = FALSE
  ) +
  geom_label(
    data = endpoint_summary,
    aes(label_x, label_y, label = label, colour = endpoint),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 5.6 / ggplot2::.pt,
    lineheight = 0.95,
    family = "Arial",
    fill = alpha("white", 0.92),
    linewidth = 0.25,
    label.padding = unit(0.15, "lines"),
    show.legend = FALSE
  ) +
  scale_x_log10(
    limits = c(0.0027, 0.27),
    breaks = c(0.003, 0.01, 0.03, 0.1),
    labels = scales::label_number(accuracy = 0.001)
  ) +
  scale_y_log10(
    limits = c(0.0027, 0.95),
    breaks = c(0.003, 0.01, 0.03, 0.1, 0.3, 0.9),
    labels = scales::label_number(accuracy = 0.001)
  ) +
  scale_colour_manual(values = endpoint_colours[c("TF–gene", "Region–gene")]) +
  scale_fill_manual(values = endpoint_colours[c("TF–gene", "Region–gene")]) +
  labs(
    title = "AP versus prevalence",
    tag = "a",
    x = "Positive-edge prevalence",
    y = "Average precision"
  ) +
  theme_genome(6.5) +
  panel_title_theme +
  theme(
    panel.grid.major.x = element_line(colour = GRID_COLOUR, linewidth = 0.25),
    plot.margin = margin(5, 5, 3, 4)
  )

auroc_points <- copy(static_region)
auroc_points[,
  x_position := as.numeric(endpoint) + seq(-0.12, 0.12, length.out = .N),
  by = endpoint
]
endpoint_summary[, x_position := as.numeric(endpoint)]

figs6b <- ggplot() +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.4
  ) +
  geom_point(
    data = auroc_points,
    aes(x_position, auroc, colour = endpoint),
    size = 1.55,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  geom_linerange(
    data = endpoint_summary,
    aes(x_position, ymin = auroc_q25, ymax = auroc_q75, colour = endpoint),
    linewidth = 1.4,
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_point(
    data = endpoint_summary,
    aes(x_position, auroc, fill = endpoint),
    shape = 21,
    size = 3.1,
    stroke = 0.75,
    colour = INK_COLOUR,
    show.legend = FALSE
  ) +
  geom_text(
    data = endpoint_summary,
    aes(x_position, auroc_q75 + 0.035, label = sprintf("%.3f", auroc)),
    size = 5.7 / ggplot2::.pt,
    family = "Arial",
    fontface = "bold",
    colour = INK_COLOUR
  ) +
  scale_x_continuous(
    limits = c(0.55, 2.45),
    breaks = 1:2,
    labels = levels(static_region$endpoint)
  ) +
  scale_y_continuous(
    limits = c(0.48, 0.88),
    breaks = seq(0.5, 0.8, by = 0.1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_colour_manual(values = endpoint_colours[c("TF–gene", "Region–gene")]) +
  scale_fill_manual(values = endpoint_colours[c("TF–gene", "Region–gene")]) +
  labs(title = "AUROC", tag = "b", x = NULL, y = "AUROC") +
  theme_genome(6.5) +
  panel_title_theme +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    plot.margin = margin(5, 4, 3, 4)
  )

truth_data <- x[
  scope %chin% scope_contract[c("tf_static", "tf_dynamic")],
  .(
    seed,
    truth_definition = fifelse(
      scope == scope_contract[["tf_static"]],
      "Static input GRN",
      "Dynamic-union truth"
    ),
    average_precision
  )
]
truth_data[,
  truth_definition := factor(
    truth_definition,
    levels = c("Static input GRN", "Dynamic-union truth")
  )
]
truth_data[, x_position := as.numeric(truth_definition)]

truth_summary <- truth_data[,
  .(
    average_precision = median(average_precision),
    q25 = quantile(average_precision, 0.25, names = FALSE),
    q75 = quantile(average_precision, 0.75, names = FALSE)
  ),
  by = truth_definition
]
truth_summary[, x_position := as.numeric(truth_definition)]

static_ap <- truth_data[
  truth_definition == "Static input GRN"
][order(seed), average_precision]
dynamic_ap <- truth_data[
  truth_definition == "Dynamic-union truth"
][order(seed), average_precision]
median_delta <- median(static_ap - dynamic_ap)

figs6c <- ggplot(truth_data, aes(x_position, average_precision, group = seed)) +
  geom_line(colour = "#B0BEC5", linewidth = 0.45, alpha = 0.85) +
  geom_point(
    aes(colour = truth_definition),
    size = 1.6,
    alpha = 0.85,
    show.legend = FALSE
  ) +
  geom_linerange(
    data = truth_summary,
    aes(x_position, ymin = q25, ymax = q75, colour = truth_definition),
    inherit.aes = FALSE,
    linewidth = 1.4,
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_point(
    data = truth_summary,
    aes(x_position, average_precision, fill = truth_definition),
    inherit.aes = FALSE,
    shape = 21,
    size = 3.1,
    stroke = 0.75,
    colour = INK_COLOUR,
    show.legend = FALSE
  ) +
  annotate(
    "label",
    x = 1.5,
    y = 0.705,
    label = sprintf("Median paired\ndifference = %.3f", median_delta),
    size = 5.6 / ggplot2::.pt,
    family = "Arial",
    colour = INK_COLOUR,
    fill = alpha("white", 0.9),
    linewidth = 0.2,
    label.padding = unit(0.12, "lines")
  ) +
  scale_x_continuous(
    limits = c(0.7, 2.3),
    breaks = 1:2,
    labels = c("Static input\nGRN", "Dynamic-union\ntruth")
  ) +
  scale_y_continuous(
    limits = c(0.27, 0.74),
    breaks = seq(0.3, 0.7, by = 0.1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_colour_manual(
    values = c(
      `Static input GRN` = endpoint_colours[["TF–gene"]],
      `Dynamic-union truth` = endpoint_colours[["Dynamic-union truth"]]
    )
  ) +
  scale_fill_manual(
    values = c(
      `Static input GRN` = endpoint_colours[["TF–gene"]],
      `Dynamic-union truth` = endpoint_colours[["Dynamic-union truth"]]
    )
  ) +
  labs(
    title = "TF–gene truth sensitivity",
    tag = "c",
    x = NULL,
    y = "Average precision"
  ) +
  theme_genome(6.5) +
  panel_title_theme +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    plot.margin = margin(5, 4, 3, 4)
  )

component <- fread(component_path)
component_methods <- c(
  full = "full_add_delete_swap__beta_rank",
  forward = "forward_only__beta_rank",
  deletion = "full_add_delete_swap__deletion_delta_bic_rank"
)

full_fit <- component[
  method == component_methods[["full"]],
  .(seed, full_auroc = auroc, full_ap = average_precision)
]
forward_fit <- component[
  method == component_methods[["forward"]],
  .(seed, forward_auroc = auroc, forward_ap = average_precision)
]
deletion_rank <- component[
  method == component_methods[["deletion"]],
  .(seed, deletion_auroc = auroc, deletion_ap = average_precision)
]
ablation <- merge(
  merge(full_fit, forward_fit, by = "seed"),
  deletion_rank,
  by = "seed"
)

ablation_colours <- c(
  `Forward-only − full` = "#7569A3",
  `Deletion-ΔBIC − coefficient` = "#5B8C32"
)

figs6d <- ggplot(ablation, aes(full_auroc, forward_auroc)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "grey55",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  geom_point(
    size = 2.0,
    shape = 21,
    stroke = 0.45,
    colour = INK_COLOUR,
    fill = ablation_colours[["Forward-only − full"]]
  ) +
  scale_x_continuous(breaks = c(0.60, 0.70, 0.80)) +
  scale_y_continuous(breaks = c(0.60, 0.70, 0.80)) +
  coord_equal(xlim = c(0.6, 0.85), ylim = c(0.6, 0.85), expand = FALSE) +
  labs(
    title = "Search ablation AUROC",
    tag = "d",
    x = "Full search (add/del/swap)\ncoefficient rank",
    y = "Forward-only\ncoefficient rank"
  ) +
  theme_genome(6.5) +
  panel_title_theme +
  theme(
    panel.grid.major.x = element_line(colour = GRID_COLOUR, linewidth = 0.25),
    plot.margin = margin(4, 4, 2, 4)
  )

figs6e <- ggplot(ablation, aes(full_ap, forward_ap)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "grey55",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  geom_point(
    size = 2.0,
    shape = 21,
    stroke = 0.45,
    colour = INK_COLOUR,
    fill = ablation_colours[["Forward-only − full"]]
  ) +
  scale_x_continuous(breaks = seq(0.3, 0.7, by = 0.1)) +
  scale_y_continuous(breaks = seq(0.3, 0.7, by = 0.1)) +
  coord_equal(xlim = c(0.28, 0.68), ylim = c(0.28, 0.68), expand = FALSE) +
  labs(
    title = "Search ablation AP",
    tag = "e",
    x = "Full search (add/del/swap)\ncoefficient rank",
    y = "Forward-only\ncoefficient rank"
  ) +
  theme_genome(6.5) +
  panel_title_theme +
  theme(
    panel.grid.major.x = element_line(colour = GRID_COLOUR, linewidth = 0.25),
    plot.margin = margin(4, 4, 2, 4)
  )

effect_data <- rbindlist(list(
  data.table(
    seed = ablation$seed,
    metric = "AUROC",
    comparison = "Forward-only − full",
    difference = ablation$forward_auroc - ablation$full_auroc
  ),
  data.table(
    seed = ablation$seed,
    metric = "AUROC",
    comparison = "Deletion-ΔBIC − coefficient",
    difference = ablation$deletion_auroc - ablation$full_auroc
  ),
  data.table(
    seed = ablation$seed,
    metric = "AP",
    comparison = "Forward-only − full",
    difference = ablation$forward_ap - ablation$full_ap
  ),
  data.table(
    seed = ablation$seed,
    metric = "AP",
    comparison = "Deletion-ΔBIC − coefficient",
    difference = ablation$deletion_ap - ablation$full_ap
  )
))
effect_data[,
  comparison := factor(comparison, levels = names(ablation_colours))
]
effect_data[, metric_y := fifelse(metric == "AUROC", 1, 0)]
effect_data[,
  group_y := metric_y +
    fifelse(comparison == "Forward-only − full", 0.12, -0.12)
]
setorder(effect_data, metric, comparison, seed)
effect_data[,
  point_y := group_y + seq(-0.07, 0.07, length.out = .N),
  by = .(metric, comparison)
]
effect_summary <- effect_data[,
  .(difference = median(difference), group_y = unique(group_y)),
  by = .(metric, comparison)
]

figs6f <- ggplot(effect_data, aes(difference, point_y, colour = comparison)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 1.55, alpha = 0.85) +
  geom_segment(
    data = effect_summary,
    aes(
      x = difference,
      xend = difference,
      y = group_y - 0.10,
      yend = group_y + 0.10,
      colour = comparison
    ),
    inherit.aes = FALSE,
    linewidth = 1.5,
    lineend = "round"
  ) +
  scale_x_continuous(
    limits = c(-0.045, 0.045),
    breaks = c(-0.04, -0.02, 0, 0.02, 0.04),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  scale_y_continuous(
    limits = c(-0.28, 1.28),
    breaks = c(0, 1),
    labels = c("AP", "AUROC")
  ) +
  scale_colour_manual(values = ablation_colours) +
  labs(
    title = "Paired metric differences",
    tag = "f",
    x = "Paired difference\n ",
    y = NULL,
    colour = NULL
  ) +
  theme_genome(6.5) +
  panel_title_theme +
  theme(
    legend.position = "none",
    aspect.ratio = 1,
    panel.grid.major.y = element_blank(),
    plot.margin = margin(4, 4, 2, 4)
  )

legend_plot <- ggplot(effect_data, aes(difference, point_y, colour = comparison)) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = ablation_colours, name = "Ablation:") +
  theme_classic(base_family = "Arial") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 6.2, face = "bold", colour = INK_COLOUR),
    legend.text = element_text(size = 5.9, colour = INK_COLOUR),
    legend.key.size = unit(0.12, "in"),
    legend.margin = margin(t = 2, b = 0)
  )
shared_legend <- cowplot::get_legend(legend_plot)

top_row <- figs6a + figs6b + figs6c + plot_layout(widths = c(0.36, 0.30, 0.34))
bottom_row <- figs6d + figs6e + figs6f + plot_layout(widths = c(0.33, 0.33, 0.34))

figs6 <- (top_row / bottom_row / shared_legend) +
  plot_layout(heights = c(2.3, 1.9, 0.25))

pdf_path <- file.path(repo, "figures", "figs6.pdf")
save_pdf(figs6, pdf_path, 6.2, 4.7)

message(
  "Saved Supplementary Fig. S6 scMultiSim endpoint and component-ablation figure."
)
