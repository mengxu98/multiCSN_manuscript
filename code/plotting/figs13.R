#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)

# Explicitly set session theme for self-contained execution.
theme_set(theme_genome(7) + theme(text = element_text(family = BASE_FAMILY)))

# ---- Parameters ----
summary_root <- file.path(repo, "results", "gse194122", "cell_budget")
source_root <- file.path(repo, "source_data", "files")
# ---------------------------------------------------------------
formal_pdf <- file.path(repo, "figures", "figs13.pdf")

paths <- c(
  overlap = file.path(summary_root, "state_network_overlap.tsv"),
  transition = file.path(summary_root, "transition_topology_summary.tsv"),
  modules = file.path(summary_root, "focal_module_support.tsv")
)
overlap <- fread(paths[["overlap"]])
transition <- fread(paths[["transition"]])
modules <- fread(paths[["modules"]])


state_labels <- c("HSC", "MK/E", "ProE", "Ery", "Norm")
transition_labels <- c("HSC to MK/E", "MK/E to ProE", "ProE to Ery", "Ery to Norm")
endpoint_labels <- c(
  `TF-gene` = "TF-gene",
  `factorized-chain` = "Factorized chain"
)
endpoint_colours <- c(`TF-gene` = "#2F6FA5", `factorized-chain` = "#C47A2C")

overlap[, `:=`(
  state_label = factor(state_labels[state_order], levels = state_labels),
  endpoint_label = factor(
    endpoint_labels[endpoint],
    levels = unname(endpoint_labels)
  )
)]
overlap[,
  x_position := as.numeric(state_label) + seq(-0.13, 0.13, length.out = .N),
  by = .(endpoint, state_order)
]

transition[, `:=`(
  transition_label = factor(
    transition_labels[transition_order],
    levels = transition_labels
  ),
  endpoint_label = factor(
    endpoint_labels[endpoint],
    levels = unname(endpoint_labels)
  )
)]
transition[,
  x_position := as.numeric(transition_label) +
    seq(-0.13, 0.13, length.out = .N),
  by = .(endpoint, transition_order)
]

panel_theme <- theme(
  text = element_text(family = "Arial"),
  plot.tag = element_text(size = 8, face = "bold", colour = INK_COLOUR),
  strip.text = element_text(size = 6.5, face = "plain"),
  panel.grid.major.x = element_blank(),
  legend.position = "none"
)

figs13a <- ggplot(overlap, aes(as.numeric(state_label), signed_jaccard)) +
  geom_boxplot(
    aes(group = interaction(endpoint_label, state_label), colour = endpoint),
    width = 0.48,
    linewidth = 0.38,
    outlier.shape = NA,
    fill = "white"
  ) +
  geom_point(
    aes(x = x_position, colour = endpoint),
    size = 1.25,
    alpha = 0.72
  ) +
  facet_wrap(~endpoint_label, nrow = 1) +
  scale_colour_manual(values = endpoint_colours) +
  scale_x_continuous(breaks = 1:5, labels = state_labels) +
  scale_y_continuous(limits = c(0, 0.60), breaks = seq(0, 0.60, 0.15)) +
  labs(tag = "a", x = NULL, y = "Signed Jaccard with all-cell fit") +
  theme_genome(6.4) +
  panel_theme +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.margin = margin(4, 5, 2, 4)
  )

figs13b <- ggplot(
  transition,
  aes(as.numeric(transition_label), composition_similarity)
) +
  geom_boxplot(
    aes(
      group = interaction(endpoint_label, transition_label),
      colour = endpoint
    ),
    width = 0.48,
    linewidth = 0.38,
    outlier.shape = NA,
    fill = "white"
  ) +
  geom_point(
    aes(x = x_position, colour = endpoint),
    size = 1.25,
    alpha = 0.72
  ) +
  facet_wrap(~endpoint_label, nrow = 1) +
  scale_colour_manual(values = endpoint_colours) +
  scale_x_continuous(breaks = 1:4, labels = transition_labels) +
  scale_y_continuous(limits = c(0.85, 1.00), breaks = seq(0.85, 1.00, 0.05)) +
  labs(tag = "b", x = NULL, y = "Edge-turnover\ncomposition similarity") +
  theme_genome(6.4) +
  panel_theme +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.margin = margin(4, 4, 2, 5)
  )

module_order <- c(
  "KLF1 (+)",
  "MYB (+)",
  "SPI1 (+)",
  "GATA2 (+)",
  "RUNX1 (+)",
  "GATA1 (+)"
)
modules[, module_id := factor(module_id, levels = module_order)]
module_full <- modules[fit_type == "all-cell"]
chain_sample <- modules[fit_type == "budget-800", .(
  median_support = median(selected_chains / module_targets),
  q25_support = quantile(selected_chains / module_targets, 0.25, names = FALSE),
  q75_support = quantile(selected_chains / module_targets, 0.75, names = FALSE)
), by = .(module_id, state_order, state)]
chain_full <- copy(module_full)
chain_full[, target_support_fraction := selected_chains / module_targets]
figs13c <- ggplot(chain_sample, aes(state_order, median_support, group = module_id)) +
  geom_ribbon(aes(ymin = q25_support, ymax = q75_support),
              fill = MULTICSN_COLOUR, alpha = 0.16, colour = NA) +
  geom_line(aes(colour = "Subsampled budget (800 cells)"), linewidth = 0.62) +
  geom_point(aes(colour = "Subsampled budget (800 cells)"), size = 1.45) +
  geom_line(data = chain_full, aes(y = target_support_fraction,
            colour = "Full cohort (all cells)"), linewidth = 0.48, linetype = 2) +
  geom_point(data = chain_full, aes(y = target_support_fraction,
             colour = "Full cohort (all cells)"), shape = 21, fill = "white",
             stroke = 0.45, size = 1.45) +
  facet_wrap(~module_id, ncol = 3, labeller = function(labels) lapply(labels, function(x) as.list(italic_gene_labels(x)))) +
  scale_colour_manual(name = NULL, values = c(
    "Subsampled budget (800 cells)" = MULTICSN_COLOUR,
    "Full cohort (all cells)" = MUTED_COLOUR)) +
  scale_x_continuous(breaks = 1:5, labels = state_labels) +
  labs(tag = "c", x = NULL, y = "Selected chains / module targets") +
  theme_genome(6.2) + panel_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 5.5),
        strip.text = element_text(size = 6), panel.spacing = unit(0.09, "in"),
        legend.position = "bottom", legend.text = element_text(size = 5.8),
        legend.key.height = unit(0.06, "in"))
figs13 <- (figs13a | figs13b) / figs13c +
  plot_layout(heights = c(0.45, 0.55))

dir.create(source_root, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(formal_pdf), recursive = TRUE, showWarnings = FALSE)

save_pdf(figs13, formal_pdf, 5.2, 3.9)
fwrite(merge(chain_sample, chain_full[, .(module_id, state_order,
  all_cell_chain_support = target_support_fraction)],
  by = c("module_id", "state_order")),
  file.path(source_root, "figs13_panel_c_focal_chain_support.tsv"), sep = "\t")
 
fwrite(overlap, file.path(source_root, "figs13_panel_a_state_overlap.tsv"), sep = "\t")
fwrite(
  transition,
  file.path(source_root, "figs13_panel_b_transition_similarity.tsv"),
  sep = "\t"
)
target_sample <- modules[fit_type == "budget-800", .(
  median_support = median(target_support_fraction), q25_support = quantile(target_support_fraction, 0.25),
  q75_support = quantile(target_support_fraction, 0.75)
), by = .(module_id, state_order, state)]
fwrite(merge(target_sample, module_full[, .(module_id, state_order, all_cell_support = target_support_fraction)],
  by = c("module_id", "state_order")), file.path(source_root, "erythroid_focal_target_support.tsv"), sep = "\t")

message("Saved Supplementary Fig. S13 cell-budget sensitivity.")
