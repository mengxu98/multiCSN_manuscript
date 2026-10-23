#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
})

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)
source(file.path(repo, "code/functions/beeline_case_views.R"), local = TRUE)
case_view <- read_beeline_case(repo, "VSC-2000-9")
truth <- case_view$truth
edges <- case_view$edges
summary <- case_view$summary
methods <- c(
  "multiCSN", "GENIE3", "GRNBoost2", "PIDC", "PPCOR", "SCODE",
  "SCORPION", "SINCERITIES", "LEAP", "GNIPLR", "SCSGL"
)
regulators <- unique(truth$regulator)
targets <- unique(truth$target)
truth_table <- truth[, .(regulator, target, weight = as.numeric(sign))]
network_tables <- setNames(lapply(methods, function(method) {
  edges[Method == method, .(
    regulator, target,
    weight = as.numeric(within_network_signed_rank)
  )]
}), methods)
summary_lookup <- setNames(summary$AUPRC, summary$Method)
figs5 <- multiCSN::plot_network_heatmap(
  c(list("Ground truth" = truth_table), network_tables),
  regulators = regulators, targets = targets,
  show_names = TRUE, show_names_position = "all",
  rect_color = "grey92", heatmap_size_lock = TRUE,
  performance_metrics = character(), truth_cell_border_lwd = 2,
  nrow = 3,
  heatmap_title = c(
    "Ground truth",
    sprintf("%s (AUPRC %.3f)", methods, summary_lookup[methods])
  ),
  heatmap_title_color = c("black", unname(method_colours[methods])),
  heatmap_palcolor = c("#2166AC", "#F7F7F7", "#B2182B"),
  legend_name = c("Truth", rep("Weight", length(methods))),
  row_title = "Regulators"
)
save_grid_pdf(
  function() {
    grid.newpage()
    grid.draw(italicise_gene_grobs(grid::forceGrob(figs5), c(regulators, targets)))
  },
  file.path(repo, "figures", "figs5.pdf"), 13, 7
)
message("Saved Supplementary Fig. 5 VSC topology with per-network legends.")
