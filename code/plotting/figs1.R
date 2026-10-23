#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(grid)
  library(ggplot2)
})

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)
source(
  file.path(repo, "code", "functions", "beeline_plot_helpers.R"),
  local = TRUE
)
root <- file.path(repo, "results", "beeline", "metrics")
synthetic_path <- file.path(root, "synthetic__other_metrics.csv")
curated_path <- file.path(root, "curated__other_metrics.csv")

figs1 <- plot_heatmaps_combined(
  as.data.frame(fread(synthetic_path)),
  as.data.frame(fread(curated_path)),
  show_significance = FALSE
)$subdataset$heatmap
metrics <- c("AUROC", "ACC", "Precision", "Recall", "F1")

format_heatmap <- function(ht, metric) {
  ht@matrix_param$width <- unit(ncol(ht@matrix) * 0.120, "in")
  ht@matrix_param$height <- unit(nrow(ht@matrix) * 0.105, "in")
  ht@row_names_param$gp <- gpar(fontsize = 7.2, lineheight = 0.9)
  ht@column_names_param$gp <- gpar(fontsize = 7.2, lineheight = 0.9)
  ht@row_title_param$gp <- gpar(fontsize = 7.2, lineheight = 0.9)
  ht@column_title <- metric
  ht@column_title_param$gp <- gpar(fontsize = 8, fontface = "plain")
  ht@matrix_legend_param$at <- switch(metric, AUROC = c(0.2, 0.6, 1), ACC = c(0.2, 0.6, 1), F1 = c(0, 0.4, 0.8), c(0, 0.5, 1))
  ht@matrix_legend_param$direction <- "horizontal"
  ht@matrix_legend_param$legend_width <- unit(0.90, "in")
  ht@matrix_legend_param$title_position <- "topleft"
  ht@matrix_legend_param$title <- metric
  ht@matrix_legend_param$title_gp <- gpar(fontsize = 7)
  ht@matrix_legend_param$labels_gp <- gpar(fontsize = 7)
  axis_param <- ht@top_annotation@anno_list[[2]]@fun@var_env$axis_param
  axis_param$gp <- gpar(fontsize = 7)
  ht@top_annotation@anno_list[[2]]@fun@var_env$axis_param <- axis_param
  ht
}
figs1@ht_list <- Map(format_heatmap, figs1@ht_list, metrics)

ap_summary <- fread(file.path(repo, "results/beeline/ap_sensitivity/method_summary.tsv"))
ap_cases <- fread(file.path(repo, "results/beeline/ap_sensitivity/case_metrics.tsv"))
method_order <- ap_summary[Scope == "all"][order(AP), Method]
ap_all <- copy(ap_cases)[, Scope := "all"]
ap_cases[, Scope := DatasetType]
ap <- rbindlist(list(ap_all, ap_cases))
ap[, Method := factor(Method, levels = method_order)]
ap[, Scope := factor(Scope, levels = c("all", "curated", "synthetic"),
                     labels = c("All 500 cases", "200 curated cases", "300 synthetic cases"))]
means <- ap[, .(AP = mean(AP)), by = .(Scope, Method)]
fill_colours <- setNames(rep("#C4C4C4", length(method_order)), method_order)
fill_colours[["multiCSN"]] <- "#D70440"
border_colours <- setNames(rep("#787878", length(method_order)), method_order)
border_colours[["multiCSN"]] <- "#A00028"

y_faces <- ifelse(method_order == "multiCSN", "bold", "plain")
y_cols <- ifelse(method_order == "multiCSN", "#D70440", "#26343E")

ap_plot <- ggplot(ap, aes(x = AP, y = Method, fill = Method, colour = Method)) +
  geom_violin(scale = "width", trim = TRUE, width = 0.82,
              linewidth = 0.28, alpha = 0.75, orientation = "y") +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white",
               colour = "#333333", linewidth = 0.30, orientation = "y", show.legend = FALSE) +
  geom_point(data = means, shape = 23, size = 1.3, fill = "#222222",
             colour = "white", stroke = 0.35, show.legend = FALSE) +
  facet_wrap(~Scope, nrow = 1) +
  scale_fill_manual(values = fill_colours) +
  scale_colour_manual(values = border_colours) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     expand = expansion(mult = c(0.015, 0.025))) +
  labs(x = "Noninterpolated average precision (AP)", y = NULL) +
  theme_classic(base_size = 7.2, base_family = BASE_FAMILY) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold", colour = "#26343E"),
    axis.text.x = element_text(size = 7.0, colour = "#26343E"),
    axis.text.y = element_text(size = 7.2, face = y_faces, colour = y_cols),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.x = element_line(colour = "#E8ECEF", linewidth = 0.25, linetype = "dashed"),
    panel.grid.major.y = element_line(colour = "#F0F0F0", linewidth = 0.25),
    panel.spacing.x = unit(4.0, "mm"),
    plot.margin = margin(2, 3, 2, 2)
  )

figure_root <- file.path(repo, "figures")
save_pdf(ap_plot, file.path(figure_root, "figs1a.pdf"), 7.6, 2.15)
save_grid_pdf(
  function() {
    grid.newpage()
    ComplexHeatmap::draw(
      figs1,
      newpage = FALSE,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "right",
      ht_gap = unit(0.26, "in"),
      merge_legends = TRUE
    )
  },
  file.path(figure_root, "figs1b.pdf"),
  9.3,
  4.6
)
source(file.path(repo, "code/functions/pdf_assembly.R"), local = TRUE)
assemble_pdf_figure(
  stem = "figs1", page = c(183, 140),
  panels = list(pdf_panel("a", 16.5, 136, 154, label_x = 6.5), pdf_panel("b", 8, 93, 171, label_x = 6.5)),
  repo = repo, figure_root = figure_root
)
