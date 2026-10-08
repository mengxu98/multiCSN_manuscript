#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(grid)
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
assert_files(c(synthetic_path, curated_path))

figs1 <- plot_heatmaps_combined(
  as.data.frame(fread(synthetic_path)),
  as.data.frame(fread(curated_path)),
  show_significance = FALSE
)$subdataset$heatmap
metrics <- c("AUROC", "ACC", "Precision", "Recall", "F1")
stopifnot(length(figs1@ht_list) == length(metrics))

format_heatmap <- function(ht, metric) {
  ht@matrix_param$width <- unit(ncol(ht@matrix) * 0.120, "in")
  ht@matrix_param$height <- unit(nrow(ht@matrix) * 0.105, "in")
  ht@row_names_param$gp <- gpar(fontsize = 7.2, lineheight = 0.9)
  ht@column_names_param$gp <- gpar(fontsize = 7.2, lineheight = 0.9)
  ht@row_title_param$gp <- gpar(fontsize = 7.2, lineheight = 0.9)
  ht@column_title <- metric
  ht@column_title_param$gp <- gpar(fontsize = 8, fontface = "plain")
  ht@matrix_legend_param$title <- metric
  ht@matrix_legend_param$title_gp <- gpar(fontsize = 6.5)
  ht@matrix_legend_param$labels_gp <- gpar(fontsize = 6.5)
  axis_param <- ht@top_annotation@anno_list[[2]]@fun@var_env$axis_param
  axis_param$gp <- gpar(fontsize = 6.5)
  ht@top_annotation@anno_list[[2]]@fun@var_env$axis_param <- axis_param
  ht
}
figs1@ht_list <- Map(format_heatmap, figs1@ht_list, metrics)

figure_root <- file.path(repo, "figures")
save_grid_pdf(
  function() {
    grid.newpage()
    ComplexHeatmap::draw(
      figs1,
      newpage = FALSE,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      ht_gap = unit(0.16, "in"),
      merge_legends = TRUE
    )
  },
  file.path(figure_root, "figs1.pdf"),
  9.3,
  3.8
)

message("Saved Supplementary Fig. 1 as one-row five-metric BEELINE figure.")
