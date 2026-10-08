#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)

screen_path <- file.path(
  repo,
  "results",
  "gse162170",
  "loci",
  "by_class",
  "eligible_loci_by_class.tsv"
)
assert_files(screen_path)
screen <- fread(screen_path)

stopifnot(
  nrow(screen) == 220L,
  all(screen$eligible),
  !anyNA(screen$class_rank),
  !anyNA(screen$median_chain_delta_bic)
)

class_levels <- c("GluN", "IN", "IPC/RG", "OPC/Glia", "Vascular")
class_colours <- c(
  GluN = "#3B6FB6",
  IN = "#8A5BA7",
  `IPC/RG` = "#D98A2B",
  `OPC/Glia` = "#4AA684",
  Vascular = "#C4474D"
)

screen[, `:=`(
  class_label = factor(class_label, levels = class_levels),
  chains = as.numeric(chains),
  regions = as.numeric(regions),
  delta_bic = as.numeric(median_chain_delta_bic)
)]

# Clamped reproducible jitter for panel b to eliminate discrete integer overplotting
set.seed(42)
screen[, `:=`(
  regions_jit = pmin(pmax(regions + runif(.N, -0.22, 0.22), 5.0), 15.0),
  chains_jit  = pmin(pmax(chains + runif(.N, -0.38, 0.38), 20.0), 60.0)
)]

focal <- screen[target %in% c("ZNF804A", "KCNC2")]
stopifnot(nrow(focal) == 2L)

# a: within-class numerical rank against the retained deletion-BIC evidence
figs9a <- ggplot(screen, aes(class_rank, delta_bic)) +
  geom_hline(
    yintercept = 3,
    linewidth = 0.35,
    linetype = "dashed",
    colour = MUTED_COLOUR
  ) +
  annotate(
    "text",
    x = 78,
    y = 2.65,
    label = "Criterion (ΔBIC ≥ 3)",
    hjust = 1,
    vjust = 0,
    size = 1.9,
    colour = MUTED_COLOUR,
    family = BASE_FAMILY
  ) +
  geom_point(
    aes(colour = class_label, size = chains),
    alpha = 0.80,
    shape = 16
  ) +
  # Focal points highlighted with halo and outline
  geom_point(
    data = focal,
    aes(class_rank, delta_bic),
    shape = 21,
    size = 3.0,
    stroke = 0.6,
    colour = "white",
    fill = NA,
    show.legend = FALSE
  ) +
  geom_point(
    data = focal,
    aes(class_rank, delta_bic, fill = class_label),
    shape = 21,
    size = 2.4,
    stroke = 0.85,
    colour = INK_COLOUR,
    show.legend = FALSE
  ) +
  geom_label_repel(
    data = focal,
    aes(class_rank, delta_bic, label = target, colour = class_label),
    fontface = "italic",
    size = 2.1,
    family = BASE_FAMILY,
    fill = alpha("white", 0.95),
    label.padding = unit(0.14, "lines"),
    label.r = unit(0.12, "lines"),
    label.size = 0.25,
    box.padding = 0.35,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.size = 0.3,
    segment.colour = MUTED_COLOUR,
    nudge_x = c(7, 6),
    nudge_y = c(0.1, 1.2),
    seed = 42,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = class_colours, guide = "none") +
  scale_fill_manual(values = class_colours) +
  scale_size_continuous(
    name = "Chains",
    range = c(0.7, 2.3),
    breaks = c(20, 40, 60)
  ) +
  scale_x_continuous(
    breaks = seq(0, 80, by = 20),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_continuous(
    breaks = seq(3, 13, by = 2),
    limits = c(2.4, 13.5),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(x = "Within-class rank", y = "Median chain deletion-ΔBIC") +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "top",
      keyheight = unit(0.11, "in"),
      keywidth = unit(0.11, "in")
    )
  ) +
  theme_genome(6.8) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.6),
    legend.margin = margin(l = -2, r = 0, t = 0, b = 0),
    plot.margin = margin(3, 4, 3, 4)
  )

# b: the selected-region and chain windows applied to every eligible locus
figs9b <- ggplot(screen, aes(regions_jit, chains_jit)) +
  annotate(
    "rect",
    xmin = 5, xmax = 15, ymin = 20, ymax = 60,
    fill = "#F2F5F8",
    colour = "#93A2B0",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 14.8,
    y = 21.5,
    label = "Eligibility window\n(5–15 reg., 20–60 ch.)",
    hjust = 1,
    vjust = 0,
    size = 1.8,
    colour = MUTED_COLOUR,
    family = BASE_FAMILY
  ) +
  geom_point(
    aes(colour = class_label),
    size = 1.1,
    alpha = 0.82,
    shape = 16
  ) +
  # Focal points highlighted with halo and outline at exact coordinates
  geom_point(
    data = focal,
    aes(regions, chains),
    shape = 21,
    size = 3.0,
    stroke = 0.6,
    colour = "white",
    fill = NA,
    show.legend = FALSE
  ) +
  geom_point(
    data = focal,
    aes(regions, chains, fill = class_label),
    shape = 21,
    size = 2.4,
    stroke = 0.85,
    colour = INK_COLOUR,
    show.legend = FALSE
  ) +
  geom_label_repel(
    data = focal,
    aes(regions, chains, label = target, colour = class_label),
    fontface = "italic",
    size = 2.1,
    family = BASE_FAMILY,
    fill = alpha("white", 0.95),
    label.padding = unit(0.14, "lines"),
    label.r = unit(0.12, "lines"),
    label.size = 0.25,
    box.padding = 0.35,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.size = 0.3,
    segment.colour = MUTED_COLOUR,
    nudge_x = c(-2.0, 2.0),
    nudge_y = c(-4.5, 5.0),
    seed = 42,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = class_colours,
    name = "Cell class"
  ) +
  scale_fill_manual(values = class_colours) +
  scale_x_continuous(
    breaks = seq(5, 15, by = 2),
    limits = c(4.3, 15.7),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq(20, 60, by = 10),
    limits = c(17, 63),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "Selected regions", y = "Chains") +
  guides(
    colour = guide_legend(
      order = 1,
      title.position = "top",
      override.aes = list(size = 1.3),
      keyheight = unit(0.11, "in"),
      keywidth = unit(0.11, "in")
    )
  ) +
  theme_genome(6.8) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.6),
    legend.margin = margin(l = -2, r = 0, t = 0, b = 0),
    plot.margin = margin(3, 4, 3, 4)
  )

figure_root <- file.path(repo, "figures")
save_pdf(figs9a, file.path(figure_root, "figs9a.pdf"), 3.2, 2.3)
save_pdf(figs9b, file.path(figure_root, "figs9b.pdf"), 3.3, 2.3)

source(file.path(repo, "code", "functions", "pdf_assembly.R"), local = TRUE)
assemble_pdf_figure(
  stem = "figs9",
  page = c(183, 74),
  panels = list(
    pdf_panel("a", 8, 66, 81),
    pdf_panel("b", 95, 66, 83)
  ),
  repo = repo,
  figure_root = figure_root
)
message("Saved Supplementary Fig. S9 complete locus-screen panels.")
