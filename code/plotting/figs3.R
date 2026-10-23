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
path <- file.path(
  repo,
  "results",
  "beeline",
  "metrics",
  "synthetic__primary_metrics.csv"
)
figs3 <- plot_heatmaps_beeline(
  as.data.frame(fread(path)),
  dataset_type = "synthetic"
)$subdataset
if (is.null(figs3)) {
  stop("Unable to construct synthetic cell-number heatmaps.")
}
save_grid_pdf(
  function() {
    grid.newpage()
    ComplexHeatmap::draw(
      figs3$heatmap,
      newpage = FALSE,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      ht_gap = unit(0.276, "in")
    )
  },
  file.path(repo, "figures", "figs3.pdf"),
  8,
  7.5
)
message("Saved Supplementary Fig. 3 synthetic cell-number heatmaps.")
