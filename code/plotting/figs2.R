#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(grid)
})

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)
source(file.path(repo, "code", "functions", "beeline_plot_helpers.R"), local = TRUE)
path <- file.path(
  repo, "results", "beeline", "metrics",
  "curated__primary_metrics.csv"
)
assert_files(path)
figs2 <- plot_heatmaps_beeline(
  as.data.frame(fread(path)),
  dataset_type = "curated"
)$dropout
if (is.null(figs2)) stop("Unable to construct curated dropout heatmaps.")
save_grid_pdf(
  function() {
    grid.newpage()
    ComplexHeatmap::draw(
      figs2$heatmap,
      newpage = FALSE, heatmap_legend_side = "right",
      annotation_legend_side = "right", ht_gap = unit(0.276, "in")
    )
  },
  file.path(repo, "figures", "figs2.pdf"),
  7, 3.6
)
message("Saved Supplementary Fig. 2 curated-dropout heatmaps.")
