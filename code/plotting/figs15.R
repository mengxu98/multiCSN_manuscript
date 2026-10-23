suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
source("code/functions/common.R")
source("code/functions/erythroid_case_plots.R")
repo <- normalizePath(getwd(), mustWork = TRUE)
retention <- fread("results/gse194122/cases/HMGA1_target_retention.tsv")
retention <- retention[set %chin% as.character(3101:3110)]
retention <- melt(retention, id.vars = "set", measure.vars = c("baseline_recall", "budget_containment"), variable.name = "measure", value.name = "fraction")
retention[, measure := factor(measure, levels = c("baseline_recall", "budget_containment"), labels = c("Reference recall", "Budget containment"))]
a <- ggplot(retention, aes(measure, fraction)) +
  geom_boxplot(width = 0.4, fill = "#EDF1F0", colour = INK_COLOUR, linewidth = 0.35, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.07, height = 0, seed = 31), colour = "#3676A3", size = 1.5) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Shared genes / set size", title = italic_gene_labels("HMGA1: smaller, largely nested loss sets")[[1L]], tag = "a") +
  theme_genome(7) + theme(plot.title = element_text(size = 7, hjust = 0.5), axis.text.x = element_text(size = 6.2))
links <- fread("results/gse194122/cases/SOX6_signed_link_recurrence.tsv")
gene_order <- sort(unique(links$target))
links[, target := factor(target, levels = rev(gene_order))]
b <- ggplot(links, aes(state_order, target, fill = budget_draws)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_tile(data = links[in_reference == TRUE], fill = NA, colour = INK_COLOUR, linewidth = 0.35) +
  geom_text(aes(label = budget_draws), size = 6 / ggplot2::.pt, family = BASE_FAMILY) +
  scale_fill_gradient(low = "#F3F5F3", high = "#3676A3", limits = c(0, 10), breaks = c(0, 5, 10), name = "Draws") +
  scale_x_continuous(breaks = 1:5, labels = c("HSC", "MK/E", "ProE", "Ery", "Norm"), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = italic_gene_labels("SOX6: positive direct links")[[1L]], tag = "b") +
  theme_genome(7) + theme(plot.title = element_text(size = 7, hjust = 0.5), axis.text.y = element_text(size = 6, face = "italic"),
    axis.text.x = element_text(size = 6), legend.title = element_text(size = 6), legend.text = element_text(size = 6),
    legend.key.height = unit(0.12, "in"), panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())
c <- erythroid_program_plot(repo, "KLF1_heme", italic_gene_labels("KLF1 / heme-associated targets (n = 13)")[[1L]]) + labs(tag = "c")
d <- erythroid_program_plot(repo, "BNC2_chain_immune", italic_gene_labels("BNC2 chain targets / immune annotation (n = 39)")[[1L]], show_tf = FALSE) + labs(tag = "d") + theme(legend.position = "none")
figure <- ((a | b) / (c | d) + plot_layout(heights = c(1.25, 1))) & theme(plot.tag = element_text(face = "bold", size = 9))
save_pdf(figure, "figures/figs15.pdf", 7.2, 5.4)
