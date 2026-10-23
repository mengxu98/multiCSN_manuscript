#!/usr/bin/env Rscript

if (!l10n_info()[["UTF-8"]]) {
  for (utf8_locale in c("C.UTF-8", "en_US.UTF-8", "UTF-8")) {
    suppressWarnings(Sys.setlocale("LC_CTYPE", utf8_locale))
    if (l10n_info()[["UTF-8"]]) break
  }
  if (!l10n_info()[["UTF-8"]]) stop("Figure 5 rendering requires a UTF-8 character locale.")
}

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code/functions/erythroid_case_plots.R"))
# Bind annotation resources to the submission snapshot. A host-wide 'latest'
# cache otherwise changes GO results while leaving gene clustering unchanged.
go_resource_root <- file.path(repo, "data", "resources", "fig5_go_3.22.0")
R.cache::setCacheRootPath(file.path(go_resource_root, "cache"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(Matrix)
  library(Seurat)
  library(HDF5Array)
  library(ComplexHeatmap)
  library(grid)
})

source(file.path(repo, "code", "functions", "common.R"), local = TRUE)
source(file.path(repo, "code", "functions", "focal_network_plot.R"), local = TRUE)
analysis_root <- file.path(
  repo,
  "results", "gse194122", "programs"
)
paths <- c(
  embedding = file.path(
    repo,
    "results/gse194122/trajectory/erythroid_gex_umap.tsv.gz"
  ),
  trajectory = file.path(
    repo,
    "results/gse194122/trajectory/erythroid_coordinates_and_pseudotime.tsv.gz"
  ),
  states = file.path(analysis_root, "module_state_summary.tsv.gz"),
  heroes = file.path(analysis_root, "hero_regulons.tsv"),
  network = file.path(analysis_root, "hero_regulon_network.tsv.gz"),
  variance = file.path(analysis_root, "module_model_variance.tsv.gz"),
  transition_topology = file.path(
    repo,
    "results/gse194122/trajectory/GSE5d_adjacent_transition_topology.tsv"
  ),
  rewiring_go = file.path(
    repo,
    "results/gse194122/donor_stability/stable_rewiring_GO_BP.tsv.gz"
  )
)
figure_root <- file.path(repo, "figures")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
states <- c("HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast")
state_labels <- c("HSC", "MK/E", "ProE", "Ery", "Norm")
state_palette <- c(
  HSC = "#3569A8",
  `MK/E prog` = "#55A07B",
  Proerythroblast = "#E4A33C",
  Erythroblast = "#D6675D",
  Normoblast = "#7966A8"
)
theme_set(theme_genome(7) + theme(text = element_text(family = BASE_FAMILY)))

if ("--panel-d" %in% commandArgs(trailingOnly = TRUE)) {
  erythroid_dynamic_panel(repo, fread(paths[["heroes"]])[order(hero_rank)],
    fread(paths[["network"]]), fread(paths[["trajectory"]])[, .(cell_id, GEX_pseudotime)],
    fread(paths[["embedding"]]), state_palette)
  source(file.path(repo, "code/functions/pdf_assembly.R"), local = TRUE)
  assemble_erythroid_main_figure(repo, figure_root)
  quit(save = "no")
}

# Compensate for cropped-panel scaling in the assembled figure.
# Final axis strokes are 0.30 mm; tick lengths are 1 mm.
panel_axis_theme <- function(display_scale) {
  theme(
    axis.line = element_line(colour = INK_COLOUR, linewidth = 0.30 / (0.75 * display_scale)),
    axis.ticks = element_line(colour = INK_COLOUR, linewidth = 0.30 / (0.75 * display_scale)),
    axis.ticks.length = unit(1 / display_scale, "mm")
  )
}

embedding <- fread(paths[["embedding"]])
trajectory <- fread(paths[["trajectory"]])[, .(cell_id, GEX_pseudotime)]
embedding <- merge(
  embedding,
  trajectory,
  by = "cell_id",
  all.x = TRUE,
  sort = FALSE
)
embedding[, state := factor(state, levels = states)]
pseudo_counts <- sparseMatrix(
  i = rep(1L, nrow(embedding)),
  j = seq_len(nrow(embedding)),
  x = 1,
  dims = c(1L, nrow(embedding)),
  dimnames = list("placeholder", embedding$cell_id)
)
umap_metadata <- as.data.frame(embedding[, .(
  state,
  donor_id,
  site,
  batch,
  sample,
  GEX_pseudotime
)])
rownames(umap_metadata) <- embedding$cell_id
umap_object <- CreateSeuratObject(pseudo_counts, meta.data = umap_metadata)
umap_matrix <- as.matrix(embedding[, .(UMAP_1, UMAP_2)])
rownames(umap_matrix) <- embedding$cell_id
colnames(umap_matrix) <- c("ERUMAP_1", "ERUMAP_2")
umap_object[["erythroid_umap"]] <- CreateDimReducObject(
  embeddings = umap_matrix,
  key = "ERUMAP_",
  assay = DefaultAssay(umap_object)
)
fig5a_result <- scop::CellDimPlot(
  umap_object,
  group.by = "state",
  reduction = "erythroid_umap",
  lineages = "GEX_pseudotime",
  lineages_trim = c(0.02, 0.98),
  lineages_span = 0.45,
  lineages_palcolor = "#6F8FAF",
  lineages_arrow = grid::arrow(length = grid::unit(0.08, "inches"), type = "closed"),
  lineages_linewidth = 1.55,
  lineages_line_bg = "white",
  lineages_line_bg_stroke = 0.65,
  label.by = "state",
  label = TRUE,
  label_insitu = TRUE,
  label_repel = TRUE,
  label.size = 4.4,
  label_repulsion = 30,
  label.fg = "#595959",
  label.bg = "white",
  label.bg.r = 0.13,
  label_point_size = 0,
  label_segment_color = "grey65",
  palcolor = unname(state_palette),
  pt.size = 0.58,
  pt.alpha = 0.55,
  raster = FALSE,
  show_stat = FALSE,
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  legend.position = "none",
  legend.direction = "horizontal",
  legend.title = "Ordered state",
  theme_use = "theme_scop",
  combine = FALSE,
  verbose = FALSE
)
fig5a <- regularise_grob_fonts(patchwork::patchworkGrob(fig5a_result[[1]]))

heroes <- fread(paths[["heroes"]])[order(hero_rank)]
network <- fread(paths[["network"]])
setorder(
  network,
  regulator,
  -represented_states,
  -stable_chains,
  -exact_regions,
  target
)
network[, display_rank := seq_len(.N), by = regulator]
shared_targets <- network[, .(regulators = uniqueN(regulator)), by = target][
  regulators > 1L,
  target
]
display_edges <- network[display_rank <= 6L | target %chin% shared_targets]
edge_df <- unique(display_edges[, .(
  from = regulator,
  to = target,
  direction,
  stable_chains
)])
nodes <- data.table(name = unique(c(edge_df$from, edge_df$to)))
nodes[
  ,
  node_type := fifelse(name %chin% heroes$regulator, "Focal TF", "Target")
]
nodes <- merge(
  nodes,
  unique(network[, .(name = target, peak_state = peak_expression_state)]),
  by = "name",
  all.x = TRUE
)
nodes <- merge(
  nodes,
  heroes[, .(
    name = regulator,
    target_genes,
    tf_peak = states[peak_gene_state_order]
  )],
  by = "name",
  all.x = TRUE
)
nodes[node_type == "Focal TF", peak_state := tf_peak]
nodes[, peak_state := factor(peak_state, levels = states)]
nodes[
  ,
  node_size := fifelse(
    node_type == "Focal TF",
    rescale(target_genes, c(5.8, 7.8)),
    NA_real_
  )
]
# Label the focal TFs plus the shared targets with the strongest multi-TF
# donor-stable support. The panel is 45 mm wide in the assembled figure, so the
# shared-target budget is capped to keep every gene name legible.
label_priority <- display_edges[, .(
  n_regulators = uniqueN(regulator),
  total_stable = sum(as.numeric(stable_chains))
), by = target]
setorder(label_priority, -n_regulators, -total_stable, target)
label_targets <- head(label_priority[n_regulators > 1L, target], 6L)
fig5b <- plot_focal_network(
  unique(display_edges[, .(
    regulator = regulator,
    target = target,
    weight = as.numeric(direction) * as.numeric(stable_chains)
  )]),
  tf_nodes = heroes$regulator,
  node_state = stats::setNames(as.character(nodes$peak_state), nodes$name),
  state_colors = state_palette,
  node_size = stats::setNames(nodes$node_size, nodes$name),
  label_nodes = c(heroes$regulator, label_targets),
  seed = 42,
  label.size = 3.5,
  edge_alpha = 0.55,
  legend_position = "none"
)

all_state_summary <- fread(paths[["states"]])
summary_state <- all_state_summary[module_id %chin% heroes$module_id]
summary_state[,
  gene_z := as.numeric(scale(median_gene_activity)),
  by = module_id
]
make_state_matrix <- function(value_col) {
  wide <- data.table::dcast(
    summary_state,
    module_id ~ state,
    value.var = value_col
  )
  mat <- as.matrix(wide[, -1L])
  rownames(mat) <- wide$module_id
  mat[rev(heroes$module_id), states, drop = FALSE]
}
gene_mat <- make_state_matrix("gene_z")
cre_mat <- make_state_matrix("regulatory_activity_z")
colnames(gene_mat) <- colnames(cre_mat) <- state_labels
plot_df <- rbind(
  data.table(
    module_id = rep(rownames(gene_mat), each = ncol(gene_mat)),
    state = rep(colnames(gene_mat), times = nrow(gene_mat)),
    z = as.vector(t(gene_mat)),
    Modality = "Gene expression (RNA)"
  ),
  data.table(
    module_id = rep(rownames(cre_mat), each = ncol(cre_mat)),
    state = rep(colnames(cre_mat), times = nrow(cre_mat)),
    z = as.vector(t(cre_mat)),
    Modality = "CRE network support (multiCSN)"
  )
)
plot_df[, state := factor(state, levels = state_labels)]
tf_order_bio <- c("RUNX1 (+)", "GATA2 (+)", "SPI1 (+)", "MYB (+)", "GATA1 (+)", "KLF1 (+)")
plot_df[, module_id := factor(module_id, levels = tf_order_bio)]
rho_df <- heroes[match(tf_order_bio, module_id), .(
  module_id,
  rho = activity_coupling_rho
)]
plot_df <- merge(plot_df, rho_df, by = "module_id")
plot_df[, panel_title := sprintf("%s  (ρ = %.2f)", module_id, rho)]
plot_df[, panel_title := factor(panel_title, levels = unique(panel_title[order(module_id)]))]

fig5c <- ggplot(plot_df, aes(x = state, y = z, group = Modality, colour = Modality)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3, linetype = "dashed") +
  geom_line(linewidth = 0.85) +
  geom_point(size = 1.7, stroke = 0.35, fill = "white", shape = 21) +
  facet_wrap(~panel_title, nrow = 2, labeller = function(labels) lapply(labels, function(x) as.list(italic_gene_labels(x)))) +
  scale_colour_manual(
    values = c(
      "Gene expression (RNA)" = "#C65445",
      "CRE network support (multiCSN)" = "#3569A8"
    ),
    name = NULL
  ) +
  scale_y_continuous(breaks = c(-2, 0, 2), limits = c(-2.2, 2.2)) +
  labs(x = NULL, y = "Within-regulon z score", title = "Expression vs. CRE-linked support") +
  theme_genome(6.8) +
  theme(
    plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 6.4, face = "bold"),
    strip.background = element_rect(fill = "#ECEFF1", colour = NA),
    axis.text.x = element_text(size = 5.8),
    legend.position = "bottom",
    legend.text = element_text(size = 6.2),
    legend.margin = margin(t = -4),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#F5F5F5", linewidth = 0.25),
    plot.margin = margin(2, 4, 2, 2)
  ) + panel_axis_theme(0.6614)

erythroid_dynamic_panel(repo, heroes, network, trajectory, embedding, state_palette)

ranked <- fread(paths[["states"]])
ranked <- data.table::dcast(
  ranked,
  module_id + regulator + direction_label ~ state,
  value.var = "regulatory_activity_z"
)
setDT(ranked)
ranked[, activity_shift := Normoblast - HSC]
ranked <- ranked[is.finite(activity_shift)]
setorder(ranked, activity_shift, module_id)
ranked[, rank := seq_len(.N)]
ranked[
  ,
  label := fifelse(
    rank == 1L | rank == .N | module_id %chin% heroes$module_id,
    module_id,
    NA_character_
  )
]
n_norm <- sum(ranked$activity_shift > 0)
n_hsc <- sum(ranked$activity_shift < 0)
shift_labels <- merge(
  ranked[!is.na(label), .(label, x = rank, y = activity_shift,
                          positive = activity_shift > 0)],
  data.table(
    label = c("ZNF382 (+)", "KLF1 (+)", "GATA1 (+)", "MYB (+)",
              "RUNX1 (+)", "GATA2 (+)", "SPI1 (+)", "ZBTB43 (+)"),
    label_x = c(65, 65, 65, 150, 95, 300, 310, 570),
    label_y = c(-2.8, -2.25, -1.45, 0.4, 1.1, 1.7, 2.3, 2.95),
    hjust = c(rep(0, 7), 1)
  ),
  by = "label", all.x = TRUE, sort = FALSE
)

fig5g <- ggplot(ranked, aes(x = rank, y = activity_shift)) +
  geom_hline(yintercept = 0, colour = "#78909C", linewidth = 0.23, linetype = "dashed") +
  geom_area(
    data = ranked[activity_shift <= 0],
    aes(x = rank, y = activity_shift),
    fill = "#3C78A8", alpha = 0.40
  ) +
  geom_area(
    data = ranked[activity_shift >= 0],
    aes(x = rank, y = activity_shift),
    fill = "#C65445", alpha = 0.40
  ) +
  geom_line(aes(colour = activity_shift > 0), linewidth = 0.36) +
  scale_colour_manual(
    values = c(`FALSE` = "#1E5F8A", `TRUE` = "#A83224"),
    guide = "none"
  ) +
  geom_point(
    data = ranked[!is.na(label)],
    aes(colour = activity_shift > 0),
    size = 0.8, shape = 19
  ) +
  geom_segment(
    data = shift_labels,
    aes(x = x, y = y, xend = label_x, yend = label_y),
    colour = "grey55", linewidth = 0.16, inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = shift_labels,
    aes(x = label_x, y = label_y, label = label, colour = positive, hjust = hjust),
    size = 5.0 / ggplot2::.pt, fontface = "italic", bg.color = "white", bg.r = 0.10,
    max.iter = 0, max.overlaps = Inf,
    box.padding = 0, point.padding = 0, segment.colour = NA,
    inherit.aes = FALSE
  ) +
  annotate(
    ggrepel::GeomTextRepel, x = 20, y = 2.7,
    label = sprintf("Norm > HSC\n(n = %d)", n_norm),
    hjust = 0, size = 5.0 / ggplot2::.pt, lineheight = 1.05,
    colour = "#A83224", bg.color = "white", bg.r = 0.10,
    max.iter = 0, max.overlaps = Inf,
    box.padding = 0, point.padding = 0, segment.colour = NA
  ) +
  annotate(
    ggrepel::GeomTextRepel, x = 320, y = -2.1,
    label = sprintf("HSC > Norm\n(n = %d)", n_hsc),
    hjust = 0, size = 5.0 / ggplot2::.pt, lineheight = 1.05,
    colour = "#1E5F8A", bg.color = "white", bg.r = 0.10,
    max.iter = 0, max.overlaps = Inf,
    box.padding = 0, point.padding = 0, segment.colour = NA
  ) +
  scale_x_continuous(
    breaks = c(0, 200, 400, 600),
    limits = c(0, 600),
    expand = expansion(mult = c(0.025, 0.055))
  ) +
  scale_y_continuous(breaks = seq(-3, 3, 1)) +
  coord_cartesian(ylim = c(-3.2, 3.2), clip = "off") +
  labs(
    x = "Regulon rank",
    y = "Support shift"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#ECEFF1", linewidth = 0.16),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6.7),
    aspect.ratio = 1,
    legend.position = "none",
    plot.margin = margin(2, 8, 3, 3)
  ) + panel_axis_theme(40 / (1.5 * 25.4))

variance <- fread(paths[["variance"]])
variance[, focal := module_id %chin% heroes$module_id]
variance <- merge(
  variance,
  heroes[, .(module_id, hero_rank)],
  by = "module_id",
  all.x = TRUE,
  sort = FALSE
)
variance[
  ,
  peak_gene_state := factor(
    states[peak_gene_state_order],
    levels = states
  )
]
state_dominant <- variance[partial_state_r2 > partial_donor_r2, .N]

fig5h <- ggplot(variance, aes(partial_donor_r2, partial_state_r2)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linewidth = 0.35,
    linetype = "22",
    colour = "#7F8C93"
  ) +
  geom_point(colour = "#AAB4BA", alpha = 0.50, size = 0.75) +
  geom_point(
    data = variance[focal == TRUE],
    aes(fill = peak_gene_state),
    shape = 21,
    colour = "white",
    stroke = 0.25,
    size = 1.75
  ) +
  geom_text_repel(
    data = variance[focal == TRUE],
    aes(label = regulator),
    fontface = "italic", size = 5.0 / ggplot2::.pt,
    bg.color = "white",
    bg.r = 0.10,
    min.segment.length = 0,
    box.padding = 0.18,
    point.padding = 0.12,
    max.overlaps = Inf,
    seed = 41,
    colour = INK_COLOUR
  ) +
  annotate(
    ggrepel::GeomTextRepel,
    x = 0.97,
    y = 0.03,
    label = sprintf("State > donor: %d/%d", state_dominant, nrow(variance)),
    hjust = 1,
    vjust = 0,
    size = 5.4 / ggplot2::.pt,
    fontface = "plain",
    colour = INK_COLOUR,
    bg.color = "white", bg.r = 0.10,
    max.iter = 0, max.overlaps = Inf,
    box.padding = 0, point.padding = 0, segment.colour = NA
  ) +
  scale_fill_manual(values = state_palette, guide = "none") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "Donor partial R²", y = "State partial R²") +
  theme(
    panel.grid = element_blank(),
    plot.margin = margin(2, 8, 3, 3)
  ) + panel_axis_theme(40 / (1.5 * 25.4))

budget_root <- file.path(repo, "results/gse194122/cell_budget")
budget_coupling <- fread(file.path(budget_root, "global_coupling_summary.tsv"))
budget_coupling[, draw := match(seed, sort(seed))]
coupling_median <- median(variance$activity_coupling_rho)

cats <- c("500 cells/state", "Full cohort")
df_range1 <- data.table(
  category = factor("500 cells/state", levels = cats),
  median_val = median(budget_coupling$median_rho),
  q25 = quantile(budget_coupling$median_rho, 0.25),
  q75 = quantile(budget_coupling$median_rho, 0.75),
  min_val = min(budget_coupling$median_rho),
  max_val = max(budget_coupling$median_rho)
)
df_pts_full <- data.table(category = factor("Full cohort", levels = cats), rho = coupling_median)
df_pts_draws <- data.table(category = factor("500 cells/state", levels = cats), rho = budget_coupling$median_rho)

set.seed(42)
fig5e <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#90A4AE", linewidth = 0.35) +
  geom_segment(aes(x = coupling_median, xend = coupling_median, y = 0.5, yend = 2.15),
               linetype = "dotted", colour = "#78909C", linewidth = 0.35) +
  geom_segment(data = df_range1,
               aes(y = category, yend = category, x = min_val, xend = max_val),
               colour = "#FFCDD2", linewidth = 3.2, lineend = "round") +
  geom_segment(data = df_range1,
               aes(y = category, yend = category, x = q25, xend = q75),
               colour = "#EF9A9A", linewidth = 3.2, lineend = "round") +
  ggbeeswarm::geom_quasirandom(data = df_pts_draws, aes(x = rho, y = category),
                               orientation = "y", width = 0.22,
                               shape = 21, fill = MULTICSN_COLOUR, colour = "white",
                               stroke = 0.35, size = 2.1) +
  geom_point(data = df_pts_full, aes(x = rho, y = category),
             shape = 23, fill = "#26343E", colour = "white",
             stroke = 0.45, size = 2.8) +
  scale_x_continuous(limits = c(-1.03, 0.22), breaks = c(-1.0, -0.71, -0.5, 0.0),
                     labels = c("-1.0", "-0.71", "-0.5", "0")) +
  scale_y_discrete(limits = cats, expand = expansion(add = c(0.48, 0.52))) +
  labs(x = expression("Median Spearman " * rho), y = NULL) +
  theme_classic(base_size = 7.0, base_family = BASE_FAMILY) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7.0, face = "bold", colour = INK_COLOUR),
    axis.text.x = element_text(size = 6.8, colour = INK_COLOUR),
    axis.title.x = element_text(size = 7.0, colour = INK_COLOUR),
    panel.grid.major.x = element_line(colour = "#F0F0F0", linewidth = 0.25),
    plot.margin = margin(3, 5, 2, 2)
  )

transition_labels <- c(
  "HSC → MK/E", "MK/E → ProE", "ProE → Ery", "Ery → Norm"
)
names(transition_labels) <- 1:4
endpoint_labels <- c(
  `factorized-chain` = "Factorized chain",
  `TF-gene` = "TF–gene"
)
topology_class_levels <- c(
  "gained", "shared_sign_stable", "shared_sign_flip", "lost"
)
transition_topology <- fread(paths[["transition_topology"]])
transition_topology[, `:=`(
  transition_label = factor(
    transition_labels[transition_order],
    levels = transition_labels
  ),
  endpoint_label = factor(
    endpoint_labels[endpoint],
    levels = endpoint_labels[c("factorized-chain", "TF-gene")]
  ),
  topology_class = factor(topology_class, levels = topology_class_levels)
)]

net_summary <- transition_topology[topology_class %in% c("gained", "lost"), .(
  net_fraction = sum(fraction * ifelse(topology_class == "gained", 1, -1)),
  gained_fraction = fraction[topology_class == "gained"],
  lost_fraction = fraction[topology_class == "lost"],
  gained_edges = edges[topology_class == "gained"],
  lost_edges = edges[topology_class == "lost"],
  union_edges = union_edges[1]
), by = .(endpoint_label, transition_label, transition_order)]

bars_data <- rbindlist(list(
  transition_topology[topology_class == "gained", .(
    endpoint_label, transition_label, transition_order,
    type = "Gained (+)",
    ymin = 0, ymax = fraction,
    fraction = fraction,
    edges = edges
  )],
  transition_topology[topology_class == "lost", .(
    endpoint_label, transition_label, transition_order,
    type = "Lost (-)",
    ymin = -fraction, ymax = 0,
    fraction = fraction,
    edges = edges
  )]
))
bars_data[, type := factor(type, levels = c("Gained (+)", "Lost (-)"))]
bars_data[type == "Gained (+)", text_y := ymax - 0.12]
bars_data[type == "Lost (-)", text_y := ymin + 0.12]

fig5f <- ggplot() +
  geom_hline(yintercept = 0, colour = "#263238", linewidth = 0.4) +
  geom_hline(yintercept = seq(-0.6, 0.6, 0.3), colour = "#ECEFF1", linewidth = 0.25) +
  geom_rect(
    data = bars_data,
    aes(xmin = as.numeric(transition_label) - 0.30,
        xmax = as.numeric(transition_label) + 0.30,
        ymin = ymin, ymax = ymax, fill = type),
    colour = "white", linewidth = 0.3
  ) +
  geom_text(
    data = bars_data[type == "Gained (+)"],
    aes(x = as.numeric(transition_label), y = text_y,
        label = sprintf("+%.0fk\n(%.0f%%)", round(edges / 1000), fraction * 100)),
    colour = "white", size = 2.0, fontface = "bold", lineheight = 0.8
  ) +
  geom_text(
    data = bars_data[type == "Lost (-)"],
    aes(x = as.numeric(transition_label), y = text_y,
        label = sprintf("-%.0fk\n(%.0f%%)", round(edges / 1000), fraction * 100)),
    colour = "white", size = 2.0, fontface = "bold", lineheight = 0.8
  ) +
  geom_line(
    data = net_summary,
    aes(x = as.numeric(transition_label), y = net_fraction, group = endpoint_label),
    colour = "#1A237E", linewidth = 0.55, linetype = "solid"
  ) +
  geom_point(
    data = net_summary,
    aes(x = as.numeric(transition_label), y = net_fraction),
    shape = 21, fill = "#3949AB", colour = "white", stroke = 0.6, size = 1.8
  ) +
  geom_label(
    data = net_summary[transition_order == 3 & endpoint_label == "Factorized chain"],
    aes(x = as.numeric(transition_label) - 0.05, y = 0.11,
        label = sprintf("+%.1f%% net", net_fraction * 100)),
    fill = "#FFFFFF", colour = "#1A237E", size = 1.9, fontface = "bold",
    label.padding = unit(0.12, "lines"), label.r = unit(0.10, "lines")
  ) +
  facet_wrap(~endpoint_label, ncol = 2) +
  scale_fill_manual(
    values = c("Gained (+)" = "#2E9F74", "Lost (-)" = "#C84E4E"),
    name = NULL
  ) +
  scale_x_continuous(breaks = 1:4, labels = transition_labels) +
  scale_y_continuous(
    limits = c(-0.55, 0.72),
    breaks = seq(-0.6, 0.6, 0.3),
    labels = function(x) paste0(ifelse(x > 0, "+", ""), scales::percent(x, accuracy = 1))
  ) +
  labs(x = NULL, y = "Turnover (% union)") +
  theme_genome(6.5) +
  theme(
    axis.text.x = element_text(angle = 22, hjust = 1, vjust = 1, size = 5.6),
    strip.text = element_text(size = 6.6, face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(t = -2, b = 0),
    legend.text = element_text(size = 5.8),
    legend.key.size = unit(0.08, "in"),
    panel.grid = element_blank(),
    plot.margin = margin(2, 4, 1, 4)
  )


case_panels <- erythroid_main_case_panels(repo)
fig5i <- case_panels$i
fig5j <- case_panels$j

aligned_gh <- cowplot::align_plots(fig5g, fig5h, align = "hv", axis = "tblr")
fig5g <- aligned_gh[[1L]]
fig5h <- aligned_gh[[2L]]

save_grid_pdf(
  function() {
    grid::grid.newpage()
    grid::grid.draw(fig5a)
  },
  file.path(figure_root, "fig5a.pdf"),
  3,
  3
)
save_pdf(fig5b, file.path(figure_root, "fig5b.pdf"), 3.3, 3)
save_pdf(fig5c, file.path(figure_root, "fig5c.pdf"), 4.15, 2.3)
save_pdf(fig5e, file.path(figure_root, "fig5e.pdf"), 2.5, 2)
save_pdf(fig5f, file.path(figure_root, "fig5f.pdf"), 4.3, 2)
save_grid_pdf(
  function() { grid::grid.newpage(); grid::grid.draw(fig5g) },
  file.path(figure_root, "fig5g.pdf"), 1.5, 1.5
)
save_grid_pdf(
  function() { grid::grid.newpage(); grid::grid.draw(fig5h) },
  file.path(figure_root, "fig5h.pdf"), 1.5, 1.5
)
save_pdf(fig5i, file.path(figure_root, "fig5i.pdf"), 3.7, 1.95)
save_pdf(fig5j, file.path(figure_root, "fig5j.pdf"), 3.7, 1.95)
message("Saved Fig. 5a-j with HMGA1 and SOX6 stage-associated programs.")

source(file.path(repo, "code", "functions", "pdf_assembly.R"), local = TRUE)
assemble_erythroid_main_figure(repo, figure_root)
