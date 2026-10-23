#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggraph)
  library(igraph)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
data.table::setDTthreads(1L)

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
result_root <- file.path(repo, "results", "gse274113", "comparison")
dgh_root <- file.path(repo, "results", "gse274113", "neighborhood")
figure_root <- file.path(repo, "figures")
analysis_root <- file.path(repo, "source_data", "files")
dir.create(analysis_root, recursive = TRUE, showWarnings = FALSE)

source(file.path(repo, "code", "functions", "common.R"), local = TRUE)
source(file.path(repo, "code", "functions", "pdf_assembly.R"), local = TRUE)

panel_a_only <- FALSE
assemble_fig3 <- function() {
  assemble_pdf_figure(
    stem = "fig3",
    page = c(183, 235),
    panels = list(
      pdf_panel("a", 8, 231, 75),
      pdf_panel("b", 8, 192, 35.5),
      pdf_panel("c", 47.5, 192, 35.5),
      pdf_panel("d", 89.2, 231, 88),
      pdf_panel("e", 8, 154, 35.5),
      pdf_panel("f", 47.5, 154, 35.5),
      pdf_panel("g", 8, 112, 81.5),
      pdf_panel("h", 96.5, 112, 78.5),
      pdf_panel("i", 8, 68, 120),
      pdf_panel("j", 130, 68, 45)
    ),
    repo = repo,
    figure_root = figure_root
  )
}

build_fig3_retention_heatmap <- function(repo_root) {
  input_path <- file.path(
    repo_root, "results", "gse274113", "stability", "layer_attrition.tsv"
  )
  layer_attrition <- data.table::fread(input_path)
  timepoints <- c("day 7", "day 9", "day 11", "day 14")
  layer_order <- c(
    "rna_detection_eligible", "rna_model_eligible",
    "multiomic_candidate_and_rna_model", "fit_eligible", "selected"
  )
  attrition_plot <- layer_attrition[layer %in% layer_order]
  attrition_plot[, `:=`(
    layer = factor(layer, levels = layer_order),
    fold_label = paste(timepoint, heldout_replicate, sep = " · "),
    replicate_number = as.integer(sub("^rep", "", heldout_replicate))
  )]
  fold_levels <- attrition_plot[
    order(match(timepoint, timepoints), replicate_number), unique(fold_label)
  ]
  attrition_wide <- data.table::dcast(
    attrition_plot,
    fold_label + timepoint + heldout_replicate + replicate_number ~ layer,
    value.var = "fraction_of_complete_edges"
  )
  attrition_wide[, fold_label := factor(fold_label, levels = fold_levels)]
  data.table::setorder(attrition_wide, fold_label)
  attrition_matrix <- as.matrix(attrition_wide[, ..layer_order])
  rownames(attrition_matrix) <- attrition_wide$heldout_replicate
  attrition_matrix <- t(attrition_matrix)
  attrition_split <- factor(attrition_wide$timepoint, levels = timepoints)
  ink <- "#26343E"
  base_family <- "Arial"
  ComplexHeatmap::Heatmap(
    attrition_matrix,
    name = "Complete-edge\nretention",
    col = circlize::colorRamp2(
      c(0, 1 / 3, 2 / 3, 1),
      c("#F5F7F8", "#C7DCE8", "#75A9C5", "#2C7FB8")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = attrition_split,
    column_gap = grid::unit(0.06, "in"),
    column_title_gp = grid::gpar(
      fontsize = 7, fontface = "plain", fontfamily = base_family
    ),
    row_names_side = "left",
    row_labels = c(
      "RNA detection", "RNA model", "Rebuilt candidate", "Fit eligible",
      "Selected"
    ),
    row_names_gp = grid::gpar(fontsize = 6.5, fontfamily = base_family),
    column_names_side = "bottom",
    column_names_centered = TRUE,
    column_names_rot = 0,
    column_names_gp = grid::gpar(fontsize = 5.8, fontfamily = base_family),
    rect_gp = grid::gpar(col = "white", lwd = 0.5),
    width = grid::unit(ncol(attrition_matrix) * 0.22, "in"),
    height = grid::unit(nrow(attrition_matrix) * 0.22, "in"),
    heatmap_legend_param = list(
      at = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0%", "25%", "50%", "75%", "100%"),
      title_gp = grid::gpar(fontsize = 7, fontface = "plain"),
      labels_gp = grid::gpar(fontsize = 6.5),
      grid_width = grid::unit(0.12, "in")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      value <- attrition_matrix[i, j]
      grid::grid.text(
        scales::percent(value, accuracy = 1), x, y,
        gp = grid::gpar(
          fontsize = 6.2,
          col = if (value >= 0.67) "white" else ink,
          fontfamily = base_family
        )
      )
    }
  )
}

build_gse274113_runx1_panels <- function(repo_root) {
  repo <- repo_root
  analysis_root <- file.path(
    repo, "results", "gse274113", "neighborhood"
  )
  paths <- c(
    hero_paths = file.path(
      analysis_root, "hero_direction_consistent_paths.tsv"
    ),
    annotations = file.path(
      repo,
      "results/gse274113/comparison/selected_chain_response_annotations.tsv.gz"
    )
  )
  theme_set(theme_genome(7) + theme(text = element_text(family = "Arial")))

  hero <- fread(paths[["hero_paths"]])
  hero_tf <- unique(hero$perturbation_regulator)
  hero_day <- unique(hero$timepoint)
  direct_edges <- hero[
    order == 1L,
    .(from = hero_tf, to = target, edge_type = "Direct", path_direction)
  ]
  root_intermediate <- unique(hero[
    order == 2L,
    .(
      from = hero_tf,
      to = selected_intermediate,
      edge_type = "Second order",
      path_direction
    )
  ])
  intermediate_target <- hero[
    order == 2L,
    .(
      from = selected_intermediate,
      to = target,
      edge_type = "Second order",
      path_direction
    )
  ]
  tree_edges <- unique(rbindlist(
    list(direct_edges, root_intermediate, intermediate_target),
    use.names = TRUE
  ))
  tree_nodes <- data.table(name = unique(c(tree_edges$from, tree_edges$to)))
  tree_nodes[
    ,
    role := fcase(
      name == hero_tf, "Perturbed TF",
      name %chin% hero$selected_intermediate, "Intermediate TF",
      name %chin% direct_edges$to, "Direct target",
      default = "Second-order target"
    )
  ]
  tree_nodes <- merge(
    tree_nodes,
    unique(hero[, .(name = target, rna_logFC)]),
    by = "name",
    all.x = TRUE
  )
  tree_graph <- graph_from_data_frame(
    tree_edges,
    directed = TRUE,
    vertices = tree_nodes
  )
  hierarchy_layout <- ggraph::create_layout(
    tree_graph,
    layout = "dendrogram",
    circular = TRUE
  )
  hierarchy_nodes <- data.table::data.table(
    name = hierarchy_layout$name,
    x = hierarchy_layout$x,
    y = hierarchy_layout$y
  )
  hierarchy_nodes <- merge(
    tree_nodes,
    hierarchy_nodes,
    by = "name",
    all.x = TRUE
  )
  hierarchy_nodes[, `Node role` := data.table::fcase(
    role == "Perturbed TF", "Perturbed TF",
    role == "Intermediate TF", "Intermediate TF",
    role == "Direct target", "Direct target",
    default = "Second-order target"
  )]
  hierarchy_nodes[, node_shape := data.table::fcase(
    role %chin% c("Perturbed TF", "Intermediate TF"), 23,
    default = 21
  )]
  hierarchy_nodes[, node_size := data.table::fcase(
    role == "Perturbed TF", 9,
    role == "Intermediate TF", 5.5,
    role == "Direct target", 3.2,
    default = 2.4
  )]
  hierarchy_nodes[, label_size := data.table::fcase(
    role == "Perturbed TF", 6.2 / ggplot2::.pt,
    role == "Intermediate TF", 5.2 / ggplot2::.pt,
    role == "Direct target", 5.0 / ggplot2::.pt,
    default = 0
  )]
  hierarchy_nodes[, label_face := data.table::fcase(
    role == "Perturbed TF", "bold.italic",
    default = "italic"
  )]
  hierarchy_edges <- data.table::copy(tree_edges)
  data.table::setnames(hierarchy_edges, "edge_type", "Path order")
  hierarchy <- thisplot::NetworkPlot(
    edge = as.data.frame(hierarchy_edges),
    node = as.data.frame(hierarchy_nodes),
    from = "from",
    to = "to",
    weight = NULL,
    node_name = "name",
    edge_group = "Path order",
    edge_palcolor = c(Direct = "#3569A8", `Second order` = "#B06B8B"),
    node_group = "Node role",
    node_palcolor = c(
      "Perturbed TF" = "#E4A33C",
      "Intermediate TF" = "#F0C987",
      "Direct target" = "#3569A8",
      "Second-order target" = "#D9E2EC"
    ),
    node_shape = "node_shape",
    node_size = "node_size",
    layout = "none",
    node_coord = c("x", "y"),
    edge_width = c(0.25, 1),
    edge_alpha = 0.5,
    label = TRUE,
    label_nodes = hierarchy_nodes[role != "Second-order target", name],
    label_size = "label_size",
    label_face = "label_face",
    label.size = 2.2,
    theme_use = "theme_void",
    legend.position = "none",
    aspect.ratio = 1
  )

  chains <- fread(paths[["annotations"]])[
    perturbation_regulator == hero_tf &
      timepoint == hero_day &
      joint_evaluable == TRUE
  ]
  chains[
    ,
    direction_consistent := gene_direction_correct &
      atac_direction_correct &
      mediated_direction_correct
  ]
  hero_targets <- unique(hero$target)
  labels <- chains[target %chin% hero_targets][
    order(rna_FDR, atac_FDR),
    .SD[1L],
    by = target
  ][1:min(12L, .N)]
  plane <- ggplot(chains, aes(atac_logFC, rna_logFC)) +
    geom_hline(yintercept = 0, colour = "grey72", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey72", linewidth = 0.3) +
    geom_point(
      data = chains[direction_consistent == FALSE],
      colour = "#AEB6BB",
      size = 0.65,
      alpha = 0.26
    ) +
    geom_point(
      data = chains[direction_consistent == TRUE],
      aes(fill = chain_direction > 0),
      shape = 21,
      colour = "white",
      stroke = 0.2,
      size = 1.5,
      alpha = 0.9
    ) +
    geom_text_repel(
      data = labels,
      aes(label = target),
      size = 5.0 / ggplot2::.pt,
      max.overlaps = Inf,
      min.segment.length = 0,
      fontface = "italic",
      seed = 20260903
    ) +
    scale_fill_manual(
      values = c(`FALSE` = "#3C78A8", `TRUE` = "#b91702"),
      labels = c("Predicted negative", "Predicted positive"),
      name = NULL,
      guide = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    labs(
      x = "Held-out CRE accessibility logFC",
      y = "Held-out target RNA logFC"
    ) +
    coord_equal() +
    theme(
      panel.grid.major = element_line(colour = GRID_COLOUR, linewidth = 0.2),
      legend.position = "bottom",
      legend.key.height = unit(5.5, "pt"),
      legend.key.spacing.y = unit(0, "pt"),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(t = -2, b = -2)
    )

  list(hierarchy = hierarchy, plane = plane)
}

if (!panel_a_only) {
  runx1_panels <- build_gse274113_runx1_panels(repo)
  fig3e <- runx1_panels$hierarchy
  fig3f <- runx1_panels$plane
  fig3i <- build_fig3_retention_heatmap(repo)
}

methods <- c("multiCSN", "SCENIC+", "Pando", "LINGER")
colours <- method_colours[methods]
shapes <- c(multiCSN = 16, `SCENIC+` = 18, Pando = 15, LINGER = 17)

theme_paper <- function(base_size = 6.4) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),
      axis.title = element_text(colour = INK_COLOUR, size = base_size),
      axis.text = element_text(colour = INK_COLOUR, size = base_size - 0.4),
      axis.line = element_line(colour = INK_COLOUR, linewidth = 0.35),
      axis.ticks = element_line(colour = INK_COLOUR, linewidth = 0.35),
      panel.grid.major.y = element_line(colour = GRID_COLOUR, linewidth = 0.25),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(colour = INK_COLOUR, face = "plain", size = base_size),
      legend.title = element_text(size = base_size - 0.3, colour = INK_COLOUR),
      legend.text = element_text(size = base_size - 0.5, colour = INK_COLOUR),
      plot.margin = margin(3, 4, 3, 4)
    )
}

metric_labels <- c(
  truth_response_coverage = "Truth-response\ncoverage",
  response_average_precision = "Raw response\nAP",
  response_ap_lift_over_prevalence = "AP lift over\nprevalence",
  response_precision_at_truth_k = "Precision at\ntruth k",
  directional_precision_at_truth_k = "Directional precision\nat truth k",
  signed_spearman_all = "Signed response\nSpearman rho"
)
reference <- fread(file.path(result_root, "stratum_metrics.tsv"))[
  method == "multiCSN" & endpoint == "TF-gene" &
    comparison_universe == "full_truth_feature_universe"
]
comparison_metrics <- fread(file.path(result_root, "stratum_metrics.tsv"))
linger <- fread(file.path(result_root, "linger_stratum_metrics.tsv"))
set.seed(20260903)
tf_effect_rows <- list()
effect <- rbindlist(lapply(c("Pando", "SCENIC+", "LINGER"), function(current) {
  other <- if (current == "LINGER") linger else comparison_metrics[
    method == current & endpoint == "TF-gene" &
      comparison_universe == "full_truth_feature_universe"
  ]
  paired <- merge(
    reference, other,
    by = c("endpoint", "comparison_universe", "regulator", "timepoint"),
    suffixes = c("_multicsn", "_comparator"), sort = TRUE
  )
  rbindlist(lapply(names(metric_labels), function(metric) {
    difference <- paired[[paste0(metric, "_comparator")]] -
      paired[[paste0(metric, "_multicsn")]]
    cluster <- data.table(regulator = paired$regulator, difference = difference)[
      is.finite(difference),
      .(tf_mean_difference = mean(difference), finite_timepoints = .N),
      by = regulator
    ]
    setorder(cluster, regulator)
    if (!nrow(cluster)) {
      return(data.table(
        metric, comparator = current, mean_difference = NA_real_,
        cluster_bootstrap_ci_lower = NA_real_, cluster_bootstrap_ci_upper = NA_real_,
        paired_strata = 0L, regulator_clusters_finite = 0L,
        wilcoxon_signed_rank_p = NA_real_
      ))
    }
    tf_effect_rows[[paste(current, metric)]] <<- cluster[, .(
      metric, comparator = current, regulator, finite_timepoints,
      difference_multicsn_minus_comparator = -tf_mean_difference
    )]
    sampled <- matrix(
      sample.int(nrow(cluster), nrow(cluster) * 5000L, replace = TRUE),
      nrow = nrow(cluster)
    )
    boot <- colMeans(matrix(cluster$tf_mean_difference[sampled], nrow = nrow(cluster)))
    data.table(
      metric, comparator = current,
      mean_difference = mean(cluster$tf_mean_difference),
      cluster_bootstrap_ci_lower = unname(quantile(boot, 0.025)),
      cluster_bootstrap_ci_upper = unname(quantile(boot, 0.975)),
      paired_strata = sum(cluster$finite_timepoints),
      regulator_clusters_finite = nrow(cluster),
      wilcoxon_signed_rank_p = suppressWarnings(
        wilcox.test(cluster$tf_mean_difference, mu = 0, exact = FALSE)$p.value
      )
    )
  }))
}))
fwrite(
  rbindlist(tf_effect_rows),
  file.path(analysis_root, "fig3a_tf_cluster_effects_source.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
effect[, holm_p := p.adjust(wilcoxon_signed_rank_p, method = "holm"), by = metric]
effect[, `:=`(
  estimate = -mean_difference,
  conf_low = -cluster_bootstrap_ci_upper,
  conf_high = -cluster_bootstrap_ci_lower,
  metric_label = factor(unname(metric_labels[metric]), levels = unname(metric_labels)),
  comparator = factor(comparator, levels = rev(c("Pando", "SCENIC+", "LINGER"))),
  holm_significant = is.finite(holm_p) & holm_p < 0.05
)]
fwrite(
  effect[, .(
    metric,
    comparator = as.character(comparator),
    paired_strata,
    regulator_clusters_finite,
    estimate,
    conf_low,
    conf_high,
    wilcoxon_signed_rank_p,
    holm_p,
    evaluable = is.finite(estimate)
  )],
  file.path(analysis_root, "fig3a_multimetric_effects_source.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fig3a <- ggplot(effect, aes(estimate, comparator, colour = comparator)) +
  geom_vline(xintercept = 0, colour = MUTED_COLOUR, linetype = "dashed", linewidth = 0.35) +
  geom_errorbar(data = effect[is.finite(estimate)],
                aes(xmin = conf_low, xmax = conf_high), width = 0,
                linewidth = 0.55, orientation = "y") +
  geom_point(data = effect[is.finite(estimate)],
             aes(shape = holm_significant), size = 1.8) +
  geom_text(
    data = effect[!is.finite(estimate)],
    aes(x = 0, y = comparator, label = "N/A", colour = comparator),
    inherit.aes = FALSE, family = "Arial", size = 2.5, hjust = -0.2
  ) +
  facet_wrap(~metric_label, scales = "free_x", ncol = 3) +
  scale_colour_manual(values = colours[c("Pando", "SCENIC+", "LINGER")], drop = FALSE) +
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 8),
    labels = c(`FALSE` = "Holm P ≥ 0.05", `TRUE` = "Holm P < 0.05"), name = NULL
  ) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01), n.breaks = 3) +
  labs(x = "multiCSN - comparator (TF-mean difference)",
       y = NULL, colour = NULL) +
  theme_paper(8.0) +
  theme(axis.text.y = element_text(size = 7.0),
        panel.spacing.x = grid::unit(0.16, "in"),
        panel.spacing.y = grid::unit(0.04, "in"), legend.position = "bottom",
        legend.box.spacing = grid::unit(4, "pt"),
        legend.margin = margin(t = 0, b = -2),
        plot.margin = margin(3, 4, 0, 4)) +
  guides(colour = "none", shape = guide_legend(override.aes = list(colour = INK_COLOUR)))

if (panel_a_only) {
  save_pdf(fig3a, file.path(figure_root, "fig3a.pdf"), 4.0, 2.0)
  assemble_fig3()
  quit(save = "no", status = 0L)
}

binding <- fread(file.path(result_root, "TF_region_binding_four_method_common_nonempty.tsv"))[
  evidence_source == "UniBind_direct_TFBS"
]
contact <- fread(file.path(result_root, "region_gene_contact_four_methods.tsv"))
joint <- fread(file.path(result_root, "chain_response_four_methods.tsv"))
fwrite(binding, file.path(analysis_root, "fig3b_tf_region_source.tsv"),
       sep = "\t", quote = FALSE, na = "NA")

binding[, method := factor(method, levels = methods)]
r_seq <- seq(0.005, 1.0, length.out = 300)
f1_levels <- c(0.10, 0.20, 0.40)
contour_dt <- rbindlist(lapply(f1_levels, function(f1) {
  p_vals <- (f1 * r_seq) / (2 * r_seq - f1)
  valid <- which(p_vals >= 0.24 & p_vals <= 0.92 & r_seq >= 0.005 & r_seq <= 1.0)
  if (length(valid)) {
    data.table(Recall = r_seq[valid], Precision = p_vals[valid], F1 = as.character(f1))
  } else NULL
}))

fig3b <- ggplot() +
  geom_line(
    data = contour_dt,
    aes(x = Recall, y = Precision, group = F1),
    colour = "grey82", linetype = "dotted", linewidth = 0.45, inherit.aes = FALSE
  ) +
  annotate("text", x = 0.19, y = 0.25, label = "F[1]==0.20", parse = TRUE,
           size = 4.4 / ggplot2::.pt, colour = "grey55", family = "Arial") +
  annotate("text", x = 0.52, y = 0.36, label = "F[1]==0.40", parse = TRUE,
           size = 4.4 / ggplot2::.pt, colour = "grey55", family = "Arial") +
  geom_point(
    data = binding,
    aes(micro_recall, micro_precision, colour = method, shape = method,
        size = predicted_regions_in_universe),
    stroke = 0.45
  ) +
  geom_label_repel(
    data = binding,
    aes(micro_recall, micro_precision, label = method, colour = method),
    family = "Arial", size = 5.6 / ggplot2::.pt, fontface = "bold",
    fill = alpha("white", 0.95),
    label.padding = unit(0.14, "lines"), label.r = unit(0.10, "lines"),
    linewidth = 0.25, min.segment.length = 0, segment.size = 0.25,
    nudge_x = c(0.15, 0.15, 0.15, 0.00),
    nudge_y = c(0.00, 0.00, 0.00, 0.13),
    box.padding = 0.20, point.padding = 0.25, seed = 42,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = colours, drop = FALSE) +
  scale_shape_manual(values = shapes, drop = FALSE) +
  scale_size_continuous(range = c(2.2, 4.5)) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.25), labels = scales::label_percent(accuracy = 1),
    limits = c(0, 1.02), expand = expansion(mult = c(0.04, 0.04))
  ) +
  scale_y_continuous(
    breaks = seq(0.3, 0.9, 0.2), labels = scales::label_percent(accuracy = 1),
    limits = c(0.22, 0.96), expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(x = "Recall", y = "Precision", colour = NULL, shape = NULL,
       size = "Predicted regions") +
  theme_paper(8.0) +
  theme(legend.position = "none")

contact[, `:=`(
  method = factor(method, levels = methods),
  timepoint = factor(timepoint, levels = c("day 7", "day 9", "day 11", "day 14"))
)]
fig3c <- ggplot(
  contact,
  aes(timepoint, enrichment_over_distance_matched_expectation,
      group = method, colour = method, shape = method)
) +
  geom_hline(yintercept = 1, colour = MUTED_COLOUR, linewidth = 0.35, linetype = "dashed") +
  geom_line(linewidth = 0.55) +
  geom_point(aes(size = target_gene_blocks), stroke = 0.4) +
  scale_colour_manual(values = colours, breaks = methods, drop = FALSE) +
  scale_shape_manual(values = shapes, breaks = methods, drop = FALSE) +
  scale_size_continuous(range = c(1.8, 4.0), breaks = c(300, 900)) +
  scale_x_discrete(expand = expansion(add = c(0.25, 0.25))) +
  scale_y_continuous(limits = c(0.55, 1.75)) +
  labs(x = NULL, y = "Observed /\ndistance-matched expected",
       colour = NULL, shape = NULL, size = "Target-gene\nblocks") +
  theme_paper(8.5) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE, order = 1,
                          override.aes = list(size = 1.8, linewidth = 0.5)),
    shape = "none",
    size = guide_legend(nrow = 1, byrow = TRUE, order = 2,
                        title.position = "left",
                        override.aes = list(size = c(1.7, 2.7)))
  ) +
  theme(
    axis.text.x = element_text(size = 7.4, angle = 0, hjust = 0.5),
    axis.title.y = element_text(size = 7.4),
    legend.position = "bottom", legend.box = "vertical",
    legend.box.just = "left", legend.justification = "left",
    legend.title = element_text(size = 7.4), legend.text = element_text(size = 7.4),
    legend.key.height = grid::unit(6, "pt"), legend.key.width = grid::unit(9, "pt"),
    legend.key.spacing.x = grid::unit(1, "pt"), legend.key.spacing.y = grid::unit(0, "pt"),
    legend.spacing.y = grid::unit(0, "pt"), legend.margin = margin(0),
    legend.box.margin = margin(0)
  )

day_levels <- c("day 7", "day 9", "day 11", "day 14")
discovery_methods <- c("multiCSN", "Pando", "SCENICplus", "LINGER")
discovery_labels <- c(
  multiCSN = "multiCSN", Pando = "Pando", SCENICplus = "SCENIC+",
  LINGER = "LINGER"
)
discovery_colours <- c(
  multiCSN = colours[["multiCSN"]], Pando = colours[["Pando"]],
  SCENICplus = colours[["SCENIC+"]], LINGER = colours[["LINGER"]]
)

formal_enrichment <- fread(file.path(
  repo, "results", "gse274113", "neighborhood", "full_background_enrichment.tsv.gz"
))[neighborhood_order == 1L]
linger_enrichment <- fread(file.path(
  dgh_root, "fig3d_linger_direct_enrichment.tsv.gz"
))
tf_order <- formal_enrichment[
  , .(median_NES = median(NES, na.rm = TRUE)),
  by = perturbation_regulator
][order(median_NES), perturbation_regulator]

direct_enrichment <- rbindlist(
  list(formal_enrichment, linger_enrichment), use.names = TRUE, fill = TRUE
)[method %chin% discovery_methods & perturbation_regulator %chin% tf_order]
direct_enrichment[, padj_global := NA_real_]
direct_enrichment[
  is.finite(pval),
  padj_global := p.adjust(pval, method = "BH"),
  by = .(modality, neighborhood_order)
]
direct_enrichment[, `:=`(
  column_id = paste(method, timepoint, sep = "|"),
  row_id = paste(modality, perturbation_regulator, sep = "|")
)]

column_ids <- unlist(lapply(discovery_methods, function(method) {
  paste(method, day_levels, sep = "|")
}))
row_ids <- c(
  paste("RNA", tf_order, sep = "|"),
  paste("ATAC", tf_order, sep = "|")
)
matrix_from_column <- function(value_column) {
  matrix_out <- matrix(
    NA_real_, nrow = length(row_ids), ncol = length(column_ids),
    dimnames = list(row_ids, column_ids)
  )
  source <- direct_enrichment[
    row_id %chin% row_ids & column_id %chin% column_ids
  ]
  index <- cbind(
    match(source$row_id, row_ids), match(source$column_id, column_ids)
  )
  matrix_out[index] <- source[[value_column]]
  matrix_out
}
nes_matrix <- matrix_from_column("NES")
padj_matrix <- matrix_from_column("padj_global")

significant_matrix <- is.finite(padj_matrix) & padj_matrix < 0.05
significant_counts <- vapply(discovery_methods, function(method) {
  columns <- startsWith(colnames(significant_matrix), paste0(method, "|"))
  rowSums(significant_matrix[, columns, drop = FALSE])
}, numeric(nrow(significant_matrix)))
colnames(significant_counts) <- unname(discovery_labels[discovery_methods])

coverage <- direct_enrichment[
  is.finite(size) & size > 0,
  .(median_target_set = median(as.numeric(size))),
  by = .(method, timepoint, modality)
]
coverage_matrix <- matrix(
  NA_real_, nrow = length(column_ids), ncol = 2L,
  dimnames = list(column_ids, c("RNA", "ATAC"))
)
coverage_index <- cbind(
  match(paste(coverage$method, coverage$timepoint, sep = "|"), column_ids),
  match(coverage$modality, c("RNA", "ATAC"))
)
coverage_matrix[coverage_index] <- log10(coverage$median_target_set)

formal_strata <- fread(file.path(
  repo, "results", "gse274113", "comparison", "stratum_metrics.tsv"
))[endpoint == "TF-gene" & comparison_universe == "full_truth_feature_universe"]
ap_source <- rbindlist(list(formal_strata, linger), use.names = TRUE, fill = TRUE)
ap_source[method == "SCENIC+", method := "SCENICplus"]
strata_ap <- ap_source[
  method %chin% discovery_methods & regulator %chin% tf_order,
  .(mean_ap = mean(response_average_precision, na.rm = TRUE)),
  by = .(method, timepoint)
]
ap_values <- unlist(lapply(discovery_methods, function(method_name) {
  strata_ap[method == method_name][match(day_levels, timepoint), mean_ap]
}), use.names = FALSE)

d_source <- merge(
  CJ(
    method = discovery_methods, perturbation_regulator = tf_order,
    timepoint = day_levels, modality = c("RNA", "ATAC"), unique = TRUE
  ),
  direct_enrichment[, .(
    method, perturbation_regulator, timepoint, modality, size, NES,
    pval, padj_global, evaluable, unavailable_reason
  )],
  by = c("method", "perturbation_regulator", "timepoint", "modality"),
  all.x = TRUE, sort = FALSE
)
fwrite(
  d_source,
  file.path(analysis_root, "fig3d_direct_enrichment_source.tsv.gz"),
  sep = "\t", quote = FALSE, na = "NA"
)

coverage_colours <- c(RNA = "#4C78A8", ATAC = "#E4A33C")
heatmap_width_in <- 1.60
heatmap_column_gap_in <- 0.10
heatmap_row_gap_in <- 0.08
heatmap_cell_in <- (
  heatmap_width_in - (length(discovery_methods) - 1L) * heatmap_column_gap_in
) / ncol(nes_matrix)
coverage_legend <- ComplexHeatmap::Legend(
  title = "Target-set size", labels = c("RNA", "ATAC"), type = "grid",
  legend_gp = grid::gpar(fill = coverage_colours, col = NA),
  grid_width = grid::unit(0.105, "in"),
  grid_height = grid::unit(0.105, "in"),
  gap = grid::unit(0.052, "in"),
  title_gap = grid::unit(0.052, "in"),
  title_gp = grid::gpar(fontsize = 5, fontfamily = "Arial"),
  labels_gp = grid::gpar(fontsize = 5, fontfamily = "Arial")
)
significance_legend <- ComplexHeatmap::Legend(
  title = "Significance", labels = "BH < 0.05",
  type = "points", pch = 16, size = grid::unit(0.04, "in"),
  legend_gp = grid::gpar(col = "#182126"),
  grid_width = grid::unit(0.105, "in"),
  grid_height = grid::unit(0.105, "in"),
  gap = grid::unit(0.052, "in"),
  title_gap = grid::unit(0.052, "in"),
  title_gp = grid::gpar(fontsize = 5, fontfamily = "Arial"),
  labels_gp = grid::gpar(fontsize = 5, fontfamily = "Arial")
)
method_legend <- ComplexHeatmap::Legend(
  title = "Method", labels = unname(discovery_labels[discovery_methods]),
  type = "grid",
  legend_gp = grid::gpar(
    fill = unname(discovery_colours[discovery_methods]), col = NA
  ),
  grid_width = grid::unit(0.105, "in"),
  grid_height = grid::unit(0.105, "in"),
  gap = grid::unit(0.052, "in"),
  title_gap = grid::unit(0.052, "in"),
  title_gp = grid::gpar(fontsize = 5, fontfamily = "Arial"),
  labels_gp = grid::gpar(fontsize = 5, fontfamily = "Arial")
)

coverage_breaks <- setdiff(pretty(range(coverage_matrix, na.rm = TRUE)), 0)
coverage_breaks <- coverage_breaks[
  coverage_breaks <= max(coverage_matrix, na.rm = TRUE)
]
ap_breaks <- setdiff(pretty(range(ap_values, na.rm = TRUE)), 0)
ap_breaks <- ap_breaks[ap_breaks <= max(ap_values, na.rm = TRUE)]
top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  "Median target-set\nsize (log10)" = ComplexHeatmap::anno_barplot(
    coverage_matrix, beside = TRUE, baseline = 0,
    ylim = c(0, max(coverage_matrix, na.rm = TRUE)), extend = 0,
    gp = grid::gpar(fill = coverage_colours, col = NA), border = FALSE,
    axis = TRUE,
    axis_param = list(
      at = coverage_breaks,
      labels = format(coverage_breaks, trim = TRUE, scientific = FALSE),
      gp = grid::gpar(fontsize = 6.8, fontfamily = "Arial")
    )
  ),
  "Mean response AP" = ComplexHeatmap::anno_barplot(
    ap_values, ylim = c(0, max(ap_values)), extend = 0,
    gp = grid::gpar(
      fill = rep(
        unname(discovery_colours[discovery_methods]), each = length(day_levels)
      ), col = NA
    ),
    border = FALSE, axis = TRUE,
    axis_param = list(
      at = ap_breaks,
      labels = format(ap_breaks, trim = TRUE, scientific = FALSE),
      gp = grid::gpar(fontsize = 6.8, fontfamily = "Arial")
    )
  ),
  annotation_height = grid::unit(c(0.45, 0.30), "in"),
  annotation_name_gp = grid::gpar(fontsize = 6.2, fontfamily = "Arial")
)
right_annotation <- ComplexHeatmap::rowAnnotation(
  "Days" = ComplexHeatmap::anno_barplot(
    significant_counts, beside = TRUE, baseline = 0, ylim = c(0, 4),
    gp = grid::gpar(
      fill = unname(discovery_colours[discovery_methods]), col = NA
    ),
    border = FALSE, axis = TRUE,
    axis_param = list(
      at = c(0, 2, 4), labels = c("0", "2", "4"),
      gp = grid::gpar(fontsize = 6.8, fontfamily = "Arial")
    )
  ),
  annotation_width = grid::unit(0.40, "in"),
  annotation_name_gp = grid::gpar(fontsize = 6.2, fontfamily = "Arial")
)

fig3d <- ComplexHeatmap::Heatmap(
  nes_matrix, name = "Direct-target NES",
  col = circlize::colorRamp2(
    c(0.6, 1.0, 1.5), c("#3C78A8", "#F7F7F7", "#b91702")
  ),
  na_col = "#E6E9EB", cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = factor(
    rep(c("RNA", "ATAC"), each = length(tf_order)), levels = c("RNA", "ATAC")
  ),
  column_split = factor(
    rep(discovery_methods, each = length(day_levels)), levels = discovery_methods,
    labels = unname(discovery_labels[discovery_methods])
  ),
  row_gap = grid::unit(heatmap_row_gap_in, "in"),
  column_gap = grid::unit(heatmap_column_gap_in, "in"),
  row_labels = rep(tf_order, 2L),
  column_labels = rep(c("d7", "d9", "d11", "d14"), length(discovery_methods)),
  row_names_gp = grid::gpar(
    fontsize = 5.8, fontface = "italic", fontfamily = "Arial"
  ),
  row_title_gp = grid::gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
  row_title_rot = 90,
  column_names_gp = grid::gpar(fontsize = 6.5, fontfamily = "Arial"),
  column_title = unname(discovery_labels[discovery_methods]),
  column_title_gp = grid::gpar(fontsize = 6.5, fontface = "plain", fontfamily = "Arial"),
  rect_gp = grid::gpar(col = "white", lwd = 0.4),
  top_annotation = top_annotation, right_annotation = right_annotation,
  width = grid::unit(heatmap_width_in, "in"),
  height = grid::unit(
    nrow(nes_matrix) * heatmap_cell_in + heatmap_row_gap_in, "in"
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (significant_matrix[i, j]) {
      grid::grid.circle(
        x, y, r = grid::unit(0.032, "in"),
        gp = grid::gpar(fill = "white", col = "white")
      )
      grid::grid.circle(
        x, y, r = grid::unit(0.020, "in"),
        gp = grid::gpar(fill = "#182126", col = "#182126")
      )
    }
  },
  heatmap_legend_param = list(
    title = "Direct-target\nNES",
    title_gp = grid::gpar(fontsize = 5, fontfamily = "Arial"),
    labels_gp = grid::gpar(fontsize = 5, fontfamily = "Arial"),
    grid_width = grid::unit(0.08, "in"),
    legend_height = grid::unit(0.62, "in"),
    title_gap = grid::unit(0.052, "in")
  )
)

formal_neighborhood <- fread(file.path(
  repo, "results", "gse274113", "neighborhood", "neighborhood_response_metrics.tsv.gz"
))
linger_neighborhood <- fread(file.path(
  dgh_root, "fig3g_linger_neighborhood_response.tsv.gz"
))
neighborhood <- rbindlist(
  list(formal_neighborhood, linger_neighborhood), use.names = TRUE, fill = TRUE
)
fwrite(
  neighborhood,
  file.path(analysis_root, "fig3g_order_response_source.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
ranked <- neighborhood[
  modality == "RNA" & neighborhood_order %in% 1:2 &
    is.finite(response_odds_ratio) & response_odds_ratio > 0,
  .(log2_odds_ratio = median(log2(response_odds_ratio))),
  by = .(method, perturbation_regulator, neighborhood_order)
]
direct_order <- ranked[
  neighborhood_order == 1L,
  .(direct = log2_odds_ratio), by = .(method, perturbation_regulator)
]
setorder(direct_order, method, direct, perturbation_regulator)
direct_order[, tf_rank := seq_len(.N), by = method]
ranked <- merge(
  ranked,
  direct_order[, .(method, perturbation_regulator, tf_rank)],
  by = c("method", "perturbation_regulator"), all.x = TRUE
)[is.finite(tf_rank)]
ranked[, `:=`(
  method = factor(method, levels = discovery_methods, labels = discovery_labels),
  order_label = factor(
    neighborhood_order, levels = 1:2, labels = c("Direct", "Second order")
  )
)]
g_labels <- ranked[
  neighborhood_order == 1L,
  unique(.SD[c(1L, .N)]), by = method
]
g_note <- data.table(
  method = factor("LINGER", levels = discovery_labels),
  label = "No distinct\nsecond-order\ntargets"
)
fig3g <- ggplot(
  ranked,
  aes(tf_rank, log2_odds_ratio, colour = order_label, group = perturbation_regulator)
) +
  geom_hline(yintercept = 0, colour = "grey72", linewidth = 0.28) +
  geom_line(colour = "#B7C0C5", linewidth = 0.34, alpha = 0.65) +
  geom_point(size = 1.35) +
  geom_text_repel(
    data = g_labels, aes(label = perturbation_regulator), colour = INK_COLOUR, fontface = "italic",
    size = 5.7 / ggplot2::.pt, min.segment.length = 0, max.overlaps = Inf,
    box.padding = 0.16, point.padding = 0.1, seed = 20260903
  ) +
  geom_text(
    data = g_note, aes(x = Inf, y = -Inf, label = label), inherit.aes = FALSE,
    family = "Arial", colour = colours[["LINGER"]],
    size = 5.7 / ggplot2::.pt, hjust = 1.05, vjust = -0.1, lineheight = 0.85
  ) +
  facet_wrap(~method, nrow = 1, scales = "free_x", drop = FALSE) +
  scale_colour_manual(
    values = c(Direct = "#263238", `Second order` = "#7E9BB4"), name = NULL
  ) +
  labs(
    x = "Perturbed TFs ranked within method",
    y = expression("Response enrichment (log"[2] * " OR)")
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_genome(7.2) +
  theme(
    text = element_text(family = "Arial"),
    panel.grid.major.y = element_line(colour = GRID_COLOUR, linewidth = 0.2),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "bottom", legend.margin = margin(t = 0, b = 0),
    legend.box.margin = margin(0), plot.margin = margin(3, 5, 1, 5)
  )

h_methods <- c("multiCSN", "Pando", "SCENIC+", "LINGER")
concordance <- rbindlist(
  list(formal_strata, linger), use.names = TRUE, fill = TRUE
)[
  endpoint == "TF-gene" & comparison_universe == "full_truth_feature_universe" &
    method %chin% h_methods,
  .(method, regulator, timepoint, signed_spearman_all)
]
concordance[, `:=`(
  method = factor(method, levels = rev(h_methods)),
  timepoint = factor(timepoint, levels = day_levels)
)]
fwrite(
  concordance[, .(
    method = as.character(method), regulator,
    timepoint = as.character(timepoint), signed_spearman_all,
    evaluable = is.finite(signed_spearman_all)
  )],
  file.path(analysis_root, "fig3h_response_concordance_source.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
finite_concordance <- concordance[is.finite(signed_spearman_all)]
concordance_summary <- finite_concordance[, .(
  n_finite = .N, n_total = uniqueN(concordance$regulator),
  median = median(signed_spearman_all),
  q1 = quantile(signed_spearman_all, 0.25, names = FALSE, type = 8),
  q3 = quantile(signed_spearman_all, 0.75, names = FALSE, type = 8)
), by = .(timepoint, method)]
fwrite(
  concordance_summary,
  file.path(analysis_root, "fig3h_response_concordance_summary.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
h_na <- data.table(
  timepoint = factor(day_levels, levels = day_levels),
  method = factor("LINGER", levels = rev(h_methods)),
  label = "N/A (unsigned)"
)
label_x <- 0.285
fig3h <- ggplot(
  finite_concordance, aes(signed_spearman_all, method, colour = method)
) +
  geom_vline(
    xintercept = 0, colour = MUTED_COLOUR, linetype = "dashed", linewidth = 0.3
  ) +
  geom_point(
    position = position_jitter(width = 0, height = 0.1, seed = 20260822),
    size = 1.1, alpha = 0.5
  ) +
  geom_errorbar(
    data = concordance_summary,
    aes(xmin = q1, xmax = q3, x = median),
    orientation = "y", width = 0, linewidth = 0.8
  ) +
  geom_point(
    data = concordance_summary, aes(x = median),
    shape = 21, fill = "white", stroke = 0.55, size = 1.8
  ) +
  geom_text(
    data = concordance_summary,
    aes(x = label_x, label = sprintf("%d/%d", n_finite, n_total)),
    hjust = 1, colour = INK_COLOUR, family = "Arial", size = 1.8
  ) +
  geom_text(
    data = h_na, aes(x = label_x, y = method, label = label),
    inherit.aes = FALSE, hjust = 1, colour = colours[["LINGER"]],
    family = "Arial", size = 1.8
  ) +
  facet_wrap(~timepoint, nrow = 2) +
  scale_colour_manual(values = colours, drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(
    limits = c(-0.1, 0.3), breaks = c(-0.1, 0, 0.2),
    labels = scales::label_number(accuracy = 0.05)
  ) +
  labs(x = "Signed response agreement (Spearman rho)", y = NULL) +
  theme_genome(6.5) +
  theme(
    text = element_text(family = "Arial"), legend.position = "none",
    panel.grid.major.x = element_line(colour = GRID_COLOUR, linewidth = 0.22),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    plot.margin = margin(2, 3, 2, 3)
  )

joint[, method := factor(method, levels = methods)]
finite_counts <- joint[, .(
  finite_tf_clusters = sum(is.finite(joint_positive_fraction)),
  total_tf_clusters = .N
), by = method]
finite_labels <- setNames(
  sprintf("%s\n%d/%d TFs", finite_counts$method,
          finite_counts$finite_tf_clusters, finite_counts$total_tf_clusters),
  finite_counts$method
)
fig3j <- ggplot(
  joint[is.finite(joint_positive_fraction)],
  aes(method, joint_positive_fraction, colour = method)
) +
  geom_violin(aes(fill = method), width = 0.8, alpha = 0.12,
              linewidth = 0.35, trim = FALSE) +
  geom_boxplot(width = 0.20, outlier.shape = NA, linewidth = 0.35, fill = "white") +
  geom_jitter(aes(size = evaluable_triplets), width = 0.10, height = 0,
              alpha = 0.75, stroke = 0.25) +
  scale_colour_manual(values = colours, drop = FALSE) +
  scale_fill_manual(values = colours, drop = FALSE) +
  scale_size_continuous(range = c(1.2, 4.2)) +
  scale_x_discrete(labels = finite_labels) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = NULL, y = "Joint ATAC + RNA\nsignificant fraction",
       colour = NULL, fill = NULL, size = "Evaluable\ntriplets") +
  theme_paper(8.0) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 6.8, lineheight = 0.95),
    axis.text.x = element_text(size = 6.6, angle = 25, hjust = 1, vjust = 1,
                                lineheight = 0.9),
    plot.margin = margin(2, 4, 2, 4)
  )

draw_fig3d <- function() {
  grid::grid.newpage()
  ComplexHeatmap::draw(
    fig3d, newpage = FALSE, heatmap_legend_side = "right",
    annotation_legend_side = "right",
    annotation_legend_list = list(
      significance_legend, method_legend, coverage_legend
    ),
    merge_legends = TRUE
  )
}

draw_fig3i <- function() {
  grid::grid.newpage()
  ComplexHeatmap::draw(fig3i, newpage = FALSE, heatmap_legend_side = "right")
}

save_pdf(fig3a, file.path(figure_root, "fig3a.pdf"), 4.0, 2.0)
save_pdf(fig3b, file.path(figure_root, "fig3b.pdf"), 1.8, 1.7)
save_pdf(fig3c, file.path(figure_root, "fig3c.pdf"), 2.0, 2.0)
save_grid_pdf(draw_fig3d, file.path(figure_root, "fig3d.pdf"), 5, 5)
save_pdf(fig3e, file.path(figure_root, "fig3e.pdf"), 2.4, 2.2)
save_pdf(fig3f, file.path(figure_root, "fig3f.pdf"), 1.9, 1.8)
save_pdf(fig3g, file.path(figure_root, "fig3g.pdf"), 3.5, 1.8)
save_pdf(fig3h, file.path(figure_root, "fig3h.pdf"), 3.0, 1.5)
save_grid_pdf(draw_fig3i, file.path(figure_root, "fig3i.pdf"), 6.55, 2.55)
save_pdf(fig3j, file.path(figure_root, "fig3j.pdf"), 2.45, 1.9)

assemble_fig3()
