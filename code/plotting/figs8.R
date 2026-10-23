suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
source("code/functions/common.R")
root <- "results/gse274113/ablation"
metrics <- fread(file.path(root, "TF_macro_metrics.tsv"))
budgets <- fread(file.path(root, "matched_budget_metrics.tsv"))
levels <- c("Direct", "Chain raw-min", "Chain percentile-min", "Chain geometric-mean", "Sign-consistent intersection")
labels <- c("Direct", "Chain\nraw min", "Chain\npercentile\nmin", "Chain\ngeometric\nmean", "Sign-\nconsistent\nintersection")
metrics[, representation := factor(representation, levels = levels)]
budgets[, representation := factor(representation, levels = levels)]
colours <- setNames(c("#2177B8", MULTICSN_COLOUR, "#8B9A9F", "#8B9A9F", "#8064A2"), levels)
paired_panel <- function(field, ylab, tag) {
  ggplot(metrics, aes(representation, .data[[field]], group = regulator)) +
    geom_line(colour = "#CCD2D5", linewidth = 0.25, alpha = 0.7) +
    geom_point(aes(colour = representation), size = 1.15, alpha = 0.85) +
    stat_summary(aes(group = representation), fun = median, geom = "point", shape = 95, size = 5, colour = INK_COLOUR) +
    scale_colour_manual(values = colours, guide = "none") +
    scale_x_discrete(labels = labels) +
    labs(x = NULL, y = ylab, tag = tag) + theme_genome(7) +
    theme(axis.text.x = element_text(size = 6), panel.grid.major.y = element_line(linewidth = 0.2))
}
a <- paired_panel("mean_ap", "Held-out RNA-response AP", "a")
b <- paired_panel("mean_response_coverage", "Response coverage", "b")
sub <- budgets[common_support_eligible == TRUE, .(precision = mean(precision), timepoints = .N), by = .(regulator, representation, k)]
counts <- unique(budgets[common_support_eligible == TRUE, .(regulator, day, k)])[, .(strata = .N, TFs = uniqueN(regulator)), by = k]
sub <- merge(sub, counts, by = "k")
sub[, budget_label := paste0("k = ", k, " (", TFs, " TFs; ", strata, " strata)")]
sub[, budget_label := factor(budget_label, levels = unique(budget_label[order(k)]))]
c <- ggplot(sub, aes(representation, precision, group = regulator)) +
  geom_line(colour = "#CCD2D5", linewidth = 0.25, alpha = 0.7) +
  geom_point(aes(colour = representation), size = 1.25) +
  stat_summary(aes(group = representation), fun = median, geom = "point", shape = 95, size = 5, colour = INK_COLOUR) +
  facet_wrap(~budget_label, nrow = 1) + scale_colour_manual(values = colours, guide = "none") +
  scale_x_discrete(labels = labels) + labs(x = NULL, y = "Matched-budget precision", tag = "c") + theme_genome(7) +
  theme(axis.text.x = element_text(size = 5.5), strip.text = element_text(size = 6.4))
figure <- ((a | b) / c + plot_layout(heights = c(1, 1))) &
  theme(plot.tag = element_text(face = "bold", size = 9))
save_pdf(figure, "figures/figs8.pdf", 7.2, 4.9)
