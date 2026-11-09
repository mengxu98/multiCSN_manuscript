suppressPackageStartupMessages({library(data.table); library(ggplot2)})
source("code/functions/common.R")
source("code/functions/pdf_assembly.R")

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
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
sensitivity <- list(
  a = paired_panel("mean_ap", "Held-out RNA-response AP", "a"),
  b = paired_panel("mean_response_coverage", "Response coverage", "b")
)
sub <- budgets[common_support_eligible == TRUE, .(precision = mean(precision), timepoints = .N), by = .(regulator, representation, k)]
counts <- unique(budgets[common_support_eligible == TRUE, .(regulator, day, k)])[, .(strata = .N, TFs = uniqueN(regulator)), by = k]
sub <- merge(sub, counts, by = "k")
sub[, budget_label := paste0("k = ", k, " (", TFs, " TFs; ", strata, " strata)")]
sub[, budget_label := factor(budget_label, levels = unique(budget_label[order(k)]))]
matched_budget <- ggplot(sub, aes(representation, precision, group = regulator)) +
  geom_line(colour = "#CCD2D5", linewidth = 0.25, alpha = 0.7) +
  geom_point(aes(colour = representation), size = 1.25) +
  stat_summary(aes(group = representation), fun = median, geom = "point", shape = 95, size = 5, colour = INK_COLOUR) +
  facet_wrap(~budget_label, nrow = 1) + scale_colour_manual(values = colours, guide = "none") +
  scale_x_discrete(labels = labels) + labs(x = NULL, y = "Matched-budget precision", tag = "c") + theme_genome(7) +
  theme(axis.text.x = element_text(size = 5.5), strip.text = element_text(size = 6.4))

assemble_perturbation_validation <- function() {
  stage <- tempfile("perturbation_validation_")
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--vanilla", "code/plotting/fig3.R", shQuote(paste0("--compact-panels=", stage))))
  stopifnot(status == 0L)

  for (id in names(sensitivity)) {
    panel <- sensitivity[[id]]
    panel$data <- metrics[is.finite(mean_ap) & is.finite(mean_response_coverage)]
    save_pdf(panel + ggplot2::labs(tag = NULL),
             file.path(stage, paste0("sensitivity_", id, ".pdf")), 3.6, 2.45)
  }
  accent <- sub[representation %in% c("Direct", "Chain geometric-mean", "Sign-consistent intersection"),
                .(precision = median(precision)), by = .(budget_label, representation)]
  panel_i <- matched_budget +
    ggplot2::geom_point(data = accent, ggplot2::aes(x = representation, y = precision, colour = representation),
                        inherit.aes = FALSE, shape = 95, size = 6.5) +
    ggplot2::labs(tag = NULL)
  save_pdf(panel_i, file.path(stage, "matched_budget.pdf"), 6.25, 1.63)

  inputs <- c(paste0("fig3", c("a", "c", "d", "g", "h", "i")), paste0("figs7", letters[1:5]))
  stopifnot(all(file.copy(file.path(repo, "figures", paste0(inputs, ".pdf")), stage)))

  assemble_pdf_figure("fig3", c(171, 197), list(
    pdf_panel("a", 8, 193, 72, label_x = 5.5, source = "fig3a.pdf"),
    pdf_panel("b", 8, 157, 30, label_x = 5.5, source = "panel_b.pdf"),
    pdf_panel("c", 42, 157, 39, label_x = 39.5, source = "panel_c.pdf"),
    pdf_panel("d", 87, 193, 76, label_x = 84.5, source = "fig3d.pdf"),
    pdf_panel("e", 8, 131.5, 35, label_x = 5.5, source = "panel_e.pdf"),
    pdf_panel("f", 46, 131.5, 34, label_x = 43.5, source = "panel_f.pdf"),
    pdf_panel("g", 8, 96, 75, label_x = 5.5, source = "fig3g.pdf"),
    pdf_panel("h", 88, 96, 76, label_x = 85.5, source = "fig3h.pdf"),
    pdf_panel("i", 8, 55, 156, label_x = 5.5, source = "matched_budget.pdf")
  ), repo, stage)

  assemble_pdf_figure("figS8", c(183, 203), list(
    pdf_panel("a", 11.5, 199, 160, source = "figs7a.pdf"),
    pdf_panel("b", 11.5, 155, 160, source = "figs7b.pdf"),
    pdf_panel("c", 11.5, 103, 48, source = "figs7c.pdf"),
    pdf_panel("d", 63, 103, 49.5, source = "figs7d.pdf"),
    pdf_panel("e", 116, 103, 55.5, source = "figs7e.pdf"),
    pdf_panel("f", 11.5, 52, 145, source = "fig3i.pdf")
  ), repo, stage)

  assemble_pdf_figure("figS9", c(171, 153), list(
    pdf_panel("a", 8, 149, 75, source = "sensitivity_a.pdf"),
    pdf_panel("b", 88, 149, 75, source = "sensitivity_b.pdf"),
    pdf_panel("c", 54, 89, 63, source = "fig3c.pdf")
  ), repo, stage)

  output <- file.path(repo, "figures", paste0(c("perturbation_validation",
                                              "perturbation_runtime_stability",
                                              "perturbation_representation_sensitivity"), ".pdf"))
  stopifnot(all(file.copy(file.path(stage, c("fig3.pdf", "figS8.pdf", "figS9.pdf")),
                         output, overwrite = TRUE)))
}

assemble_perturbation_validation()
