#!/usr/bin/env Rscript

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
# Bind annotation resources to the submission snapshot. A host-wide 'latest'
# cache otherwise changes GO results while leaving gene clustering unchanged.
go_resource_root <- file.path(repo, "data", "resources", "fig5_go_3.22.0")
stopifnot(dir.exists(go_resource_root))
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
source(file.path(repo, "code", "plotting", "focal_network_plot.R"), local = TRUE)
analysis_root <- file.path(
  repo,
  "results",
  "gse194122",
  "state_regulomes"
)
paths <- c(
  embedding = file.path(
    repo,
    "results/gse194122/erythroid_gex_umap/erythroid_gex_umap.tsv.gz"
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
    "results/gse194122/state_transitions/transition_summary/GSE5d_adjacent_transition_topology.tsv"
  ),
  rewiring_go = file.path(
    repo,
    "results/gse194122/donor_stable_rewiring/stable_rewiring_GO_BP.tsv.gz"
  ),
  h5ad = file.path(
    repo,
    "data/gse194122/processed_multiome.h5ad"
  )
)
# The expression layer is read either from the released H5AD (execution host) or from a
# compact slice of it that is written on first use and kept beside the results.
expression_subset <- file.path(repo, "data", "gse194122", "erythroid_expression_subset.rds")
if (file.exists(expression_subset)) {
  assert_files(paths[names(paths) != "h5ad"])
} else {
  assert_files(paths)
}
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
stopifnot(
  nrow(embedding) == 10952L,
  !anyDuplicated(embedding$cell_id),
  sum(is.finite(embedding$GEX_pseudotime)) == 10862L
)
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
regularise_grob_fonts <- function(grob) {
  if (!inherits(grob, "grob")) {
    return(grob)
  }
  if (!is.null(grob$gp)) {
    grob$gp[["font"]] <- 1L
  }
  if (inherits(grob, "textrepeltree") && "fontface" %in% names(grob$data)) {
    grob$data$fontface <- "plain"
  }
  if (!is.null(grob$children)) {
    grob$children <- do.call(
      grid::gList,
      lapply(grob$children, regularise_grob_fonts)
    )
  }
  if (!is.null(grob$grobs)) {
    grob$grobs <- lapply(grob$grobs, regularise_grob_fonts)
  }
  if (inherits(grob, "recordedGrob") && !is.null(grob$list)) {
    grob$list <- lapply(
      grob$list,
      function(x) {
        if (inherits(x, "grob")) regularise_grob_fonts(x) else x
      }
    )
  }
  grob
}
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
stopifnot(
  identical(dim(gene_mat), c(6L, 5L)),
  identical(dim(cre_mat), c(6L, 5L))
)
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
  facet_wrap(~panel_title, nrow = 2) +
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

dynamic_features <- unique(c(
  heroes$regulator,
  network[order(-represented_states, -stable_chains, target), target]
))
stopifnot(length(dynamic_features) == 446L)
dynamic_meta <- merge(
  trajectory[is.finite(GEX_pseudotime), .(cell_id, GEX_pseudotime)],
  embedding[, .(cell_id, state)],
  by = "cell_id",
  all.x = TRUE,
  sort = FALSE
)
stopifnot(nrow(dynamic_meta) == 10862L, !anyNA(dynamic_meta$state))
dynamic_matrix <- if (file.exists(expression_subset)) {
  cached <- readRDS(expression_subset)
  stopifnot(
    all(dynamic_features %in% rownames(cached)),
    all(dynamic_meta$cell_id %in% colnames(cached))
  )
  as(cached[dynamic_features, dynamic_meta$cell_id, drop = FALSE], "CsparseMatrix")
} else {
  h5_matrix <- HDF5Array::H5ADMatrix(paths[["h5ad"]])
  stopifnot(all(dynamic_features %in% rownames(h5_matrix)))
  extracted <- as(
    h5_matrix[
      match(dynamic_features, rownames(h5_matrix)),
      match(dynamic_meta$cell_id, colnames(h5_matrix)),
      drop = FALSE
    ],
    "CsparseMatrix"
  )
  colnames(extracted) <- dynamic_meta$cell_id
  dir.create(dirname(expression_subset), recursive = TRUE, showWarnings = FALSE)
  saveRDS(extracted, expression_subset)
  extracted
}
colnames(dynamic_matrix) <- dynamic_meta$cell_id
dynamic_metadata <- as.data.frame(dynamic_meta[, .(GEX_pseudotime, state)])
rownames(dynamic_metadata) <- dynamic_meta$cell_id
dynamic_object <- CreateSeuratObject(
  counts = dynamic_matrix,
  meta.data = dynamic_metadata
)
LayerData(dynamic_object, assay = "RNA", layer = "data") <- dynamic_matrix
ComplexHeatmap::ht_opt(
  legend_title_gp = grid::gpar(fontsize = 12.5, fontface = "plain"),
  legend_labels_gp = grid::gpar(fontsize = 11.5),
  heatmap_row_title_gp = grid::gpar(fontsize = 11.5),
  heatmap_column_title_gp = grid::gpar(fontsize = 12.5)
)
dynamic_result <- scop::DynamicHeatmap(
  object = dynamic_object,
  lineages = "GEX_pseudotime",
  features = dynamic_features,
  border = FALSE,
  heatmap_border = TRUE,
  heatmap_border_palcolor = "#9A9A9A",
  cell_annotation_border = FALSE,
  feature_annotation_border = FALSE,
  use_fitted = TRUE,
  layer = "data",
  assay = "RNA",
  lib_normalize = FALSE,
  exp_method = "zscore",
  exp_legend_title = "Z-score",
  cell_bins = 100,
  n_split = 5,
  split_method = "kmeans-peaktime",
  seed = 11,
  anno_terms = TRUE,
  species = "Homo_sapiens",
  db = "GO_BP",
  db_version = "3.22.0",
  minGSSize = 10,
  maxGSSize = 500,
  padjustCutoff = 0.05,
  topTerm = 5,
  terms_stat = "none",
  terms_text_palcolor = INK_COLOUR,
  terms_stat_label_palcolor = INK_COLOUR,
  terms_group_background = TRUE,
  terms_background_palcolor = "#F7F8F9",
  terms_border = TRUE,
  terms_border_palcolor = "#7A858B",
  terms_border_size = 0.8,
  nlabel = 0,
  features_label = heroes$regulator,
  cell_annotation = "state",
  cell_annotation_palcolor = list(state_palette),
  pseudotime_palcolor = c("#3569A8", "#F6C85F", "#C65445"),
  feature_split_palcolor = c(
    "#3569A8",
    "#55A07B",
    "#E4A33C",
    "#D6675D",
    "#7966A8"
  ),
  use_raster = TRUE,
  width = 1,
  height = 2.5
)
ComplexHeatmap::ht_opt(RESET = TRUE)
tf_cluster <- data.table(
  regulator = heroes$regulator,
  cluster = as.character(dynamic_result$feature_split[heroes$regulator]),
  hero_rank = heroes$hero_rank
)
stopifnot(!anyNA(tf_cluster$cluster))
setorder(tf_cluster, cluster, hero_rank)
dynamic_source_root <- file.path(
  repo, "source_data", "files"
)
dir.create(dynamic_source_root, recursive = TRUE, showWarnings = FALSE)
dynamic_feature_programs <- data.table(
  gene = names(dynamic_result$feature_split),
  cluster = as.character(dynamic_result$feature_split)
)
dynamic_feature_programs[
  ,
  feature_role := fifelse(
    gene %chin% heroes$regulator,
    "Focal TF",
    "Target"
  )
]
dynamic_feature_programs <- merge(
  dynamic_feature_programs,
  network[,
    .(
      represented_states = max(represented_states),
      stable_chains = sum(stable_chains)
    ),
    by = .(gene = target)
  ],
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
fwrite(
  dynamic_feature_programs,
  file.path(dynamic_source_root, "dynamic_feature_programs.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)
fwrite(
  as.data.table(dynamic_result$enrichment$enrichment),
  file.path(dynamic_source_root, "dynamic_go_enrichment.tsv.gz"),
  sep = "\t",
  quote = FALSE,
  na = "NA",
  compress = "gzip"
)
fig5d <- regularise_grob_fonts(attr(dynamic_result$plot, "grobs")$full)
dynamic_size <- attr(dynamic_result$plot, "size")
dynamic_size_numeric <- unlist(
  dynamic_size[c("width", "height")],
  use.names = FALSE
)

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
stopifnot(!anyNA(shift_labels$label_x))

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
    size = 5.0 / ggplot2::.pt, bg.color = "white", bg.r = 0.10,
    max.iter = 0, max.overlaps = Inf,
    box.padding = 0, point.padding = 0, segment.colour = NA,
    inherit.aes = FALSE
  ) +
  annotate(
    ggrepel::GeomTextRepel, x = 20, y = 2.8,
    label = sprintf("Higher in HSC\n(%d regulons)", n_hsc),
    hjust = 0, size = 5.0 / ggplot2::.pt, lineheight = 1.05,
    colour = "#1E5F8A", bg.color = "white", bg.r = 0.10,
    max.iter = 0, max.overlaps = Inf,
    box.padding = 0, point.padding = 0, segment.colour = NA
  ) +
  annotate(
    ggrepel::GeomTextRepel, x = 320, y = -2.1,
    label = sprintf("Higher in\nNormoblast\n(%d regulons)", n_norm),
    hjust = 0, size = 5.0 / ggplot2::.pt, lineheight = 1.05,
    colour = "#A83224", bg.color = "white", bg.r = 0.10,
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
stopifnot(
  nrow(variance) == 582L,
  all(is.finite(variance$partial_state_r2)),
  all(is.finite(variance$partial_donor_r2)),
  all(is.finite(variance$activity_coupling_rho))
)
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
state_median <- median(variance$partial_state_r2)
donor_median <- median(variance$partial_donor_r2)

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
    size = 5.0 / ggplot2::.pt,
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

coupling_median <- median(variance$activity_coupling_rho)
neg_pct <- mean(variance$activity_coupling_rho < 0) * 100

dens <- density(variance$activity_coupling_rho, from = -1, to = 1, n = 512)
dens_df <- data.frame(x = dens$x, y = dens$y)

fig5e <- ggplot() +
  geom_vline(xintercept = 0, colour = "#90A4AE", linetype = "dashed", linewidth = 0.35) +
  geom_area(data = dens_df[dens_df$x <= 0, ], aes(x = x, y = y), fill = "#3F51B5", alpha = 0.25) +
  geom_area(data = dens_df[dens_df$x >= 0, ], aes(x = x, y = y), fill = "#FF7043", alpha = 0.25) +
  geom_line(data = dens_df, aes(x = x, y = y), colour = "#263238", linewidth = 0.5) +
  geom_jitter(data = variance[focal == FALSE], aes(x = activity_coupling_rho, y = -0.06),
              height = 0.02, alpha = 0.30, size = 0.45, colour = "#78909C") +
  geom_point(data = variance[focal == TRUE], aes(x = activity_coupling_rho, y = -0.06),
             colour = "#D81B60", size = 1.3, shape = 17) +
  geom_vline(xintercept = coupling_median, colour = "#D81B60", linewidth = 0.6, linetype = "solid") +
  annotate("label", x = coupling_median + 0.05, y = max(dens_df$y) * 0.90,
           label = sprintf("Median ρ = %.3f\n(%.0f%% negative)", coupling_median, neg_pct),
           hjust = 0, size = 1.8, fontface = "bold", fill = "white", colour = "#D81B60",
           label.padding = unit(0.14, "lines"), label.r = unit(0.10, "lines"), lineheight = 0.8) +
  annotate("text", x = -0.75, y = max(dens_df$y) * 0.40, label = "Negative\ncorrelation",
           size = 1.7, fontface = "italic", colour = "#1A237E") +
  annotate("text", x = 0.55, y = max(dens_df$y) * 0.25, label = "Positive\ncorrelation",
           size = 1.7, fontface = "italic", colour = "#BF360C") +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), limits = c(-1, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.05))) +
  labs(x = "RNA–CRE support Spearman rho", y = "Density") +
  theme_genome(6.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(2, 6, 2, 4)
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
topology_class_labels <- c(
  lost = "Lost",
  shared_sign_flip = "Shared, sign flip",
  shared_sign_stable = "Shared, sign stable",
  gained = "Gained"
)
topology_class_colours <- c(
  lost = "#C84E4E",
  shared_sign_flip = "#D99A2B",
  shared_sign_stable = "#7897C1",
  gained = "#2E9F74"
)
transition_topology <- fread(paths[["transition_topology"]])
stopifnot(
  nrow(transition_topology) == 32L,
  uniqueN(transition_topology$endpoint) == 2L,
  uniqueN(transition_topology$transition_order) == 4L,
  all(abs(
    transition_topology[, sum(fraction), by = .(endpoint, transition_order)]$V1 - 1
  ) < 1e-10)
)
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

candidate_state <- copy(
  all_state_summary[regulator %chin% c("HMGA1", "BNC2")]
)
candidate_state[, rna_z := as.numeric(scale(median_gene_activity)), by = module_id]
candidate_long <- data.table::melt(
  candidate_state[, .(
    regulator,
    module_id,
    state_order,
    state,
    `RNA module` = rna_z,
    `CRE-linked support` = regulatory_activity_z
  )],
  id.vars = c("regulator", "module_id", "state_order", "state"),
  variable.name = "signal",
  value.name = "z_score"
)
candidate_signal_colours <- c(
  `RNA module` = "#D2644D",
  `CRE-linked support` = "#3676A3"
)
candidate_transition_colours <- c(gained = "#3B8F77", lost = "#C65445")

make_candidate_trajectory <- function(regulator_name) {
  panel_data <- candidate_long[regulator == regulator_name]
  peaks <- panel_data[, .SD[which.max(z_score)], by = signal]
  p <- ggplot(
    panel_data,
    aes(state_order, z_score, colour = signal, group = signal)
  ) +
    geom_hline(yintercept = 0, colour = "grey78", linewidth = 0.30) +
    geom_line(linewidth = 0.70) +
    geom_point(size = 1.65, stroke = 0) +
    scale_colour_manual(values = candidate_signal_colours, name = NULL) +
    scale_x_continuous(breaks = 1:5, labels = state_labels) +
    coord_cartesian(ylim = c(-1.98, 2.08), clip = "off") +
    labs(x = NULL, y = "Within-regulon z score") +
    theme_genome(7.0) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 5.7),
      axis.text.y = element_text(size = 5.7),
      axis.title.y = element_text(size = 6.1),
      legend.position = "none",
      legend.direction = "horizontal",
      legend.text = element_text(size = 5.5),
      legend.key.width = unit(0.13, "in"),
      legend.margin = margin(0),
      legend.box.margin = margin(0),
      plot.margin = margin(1, 3, 3, 2)
    )

  if (regulator_name == "HMGA1") {
    p <- p +
      annotate(
        ggrepel::GeomTextRepel,
        x = 3, y = peaks[signal == "RNA module", z_score] + 0.15,
        label = "RNA peak", colour = candidate_signal_colours[["RNA module"]],
        hjust = 0.55, size = 5.4 / ggplot2::.pt, family = BASE_FAMILY,
        bg.color = "white", bg.r = 0.10,
        max.iter = 0, max.overlaps = Inf,
        box.padding = 0, point.padding = 0, segment.colour = NA
      ) +
      annotate(
        ggrepel::GeomTextRepel,
        x = 4, y = peaks[signal == "CRE-linked support", z_score] + 0.15,
        label = "Support peak", colour = candidate_signal_colours[["CRE-linked support"]],
        hjust = 0.45, size = 5.4 / ggplot2::.pt, family = BASE_FAMILY,
        bg.color = "white", bg.r = 0.10,
        max.iter = 0, max.overlaps = Inf,
        box.padding = 0, point.padding = 0, segment.colour = NA
      ) +
      annotate(
        "segment",
        x = 4.05, xend = 4.91, y = -1.84, yend = -1.84,
        colour = candidate_transition_colours[["lost"]], linewidth = 0.58,
        arrow = arrow(length = unit(0.055, "in"), type = "closed")
      ) +
      annotate(
        ggrepel::GeomTextRepel,
        x = 4.48, y = -1.62, label = "targets lost",
        colour = candidate_transition_colours[["lost"]],
        size = 5.1 / ggplot2::.pt, family = BASE_FAMILY,
        bg.color = "white", bg.r = 0.10,
        max.iter = 0, max.overlaps = Inf,
        box.padding = 0, point.padding = 0, segment.colour = NA
      )
  } else {
    p <- p +
      annotate(
        ggrepel::GeomTextRepel,
        x = 1, y = peaks[signal == "RNA module", z_score] + 0.15,
        label = "RNA peak", colour = candidate_signal_colours[["RNA module"]],
        hjust = 0, size = 5.4 / ggplot2::.pt, family = BASE_FAMILY,
        bg.color = "white", bg.r = 0.10,
        max.iter = 0, max.overlaps = Inf,
        box.padding = 0, point.padding = 0, segment.colour = NA
      ) +
      annotate(
        ggrepel::GeomTextRepel,
        x = 4, y = peaks[signal == "CRE-linked support", z_score] + 0.15,
        label = "Support peak", colour = candidate_signal_colours[["CRE-linked support"]],
        hjust = 0.45, size = 5.4 / ggplot2::.pt, family = BASE_FAMILY,
        bg.color = "white", bg.r = 0.10,
        max.iter = 0, max.overlaps = Inf,
        box.padding = 0, point.padding = 0, segment.colour = NA
      ) +
      annotate(
        "segment",
        x = 3.05, xend = 3.91, y = -1.84, yend = -1.84,
        colour = candidate_transition_colours[["gained"]], linewidth = 0.58,
        arrow = arrow(length = unit(0.055, "in"), type = "closed")
      ) +
      annotate(
        ggrepel::GeomTextRepel,
        x = 3.48, y = -1.62, label = "gained",
        colour = candidate_transition_colours[["gained"]],
        size = 5.1 / ggplot2::.pt, family = BASE_FAMILY,
        bg.color = "white", bg.r = 0.10,
        max.iter = 0, max.overlaps = Inf,
        box.padding = 0, point.padding = 0, segment.colour = NA
      ) +
      annotate(
        "segment",
        x = 4.05, xend = 4.91, y = -1.84, yend = -1.84,
        colour = candidate_transition_colours[["lost"]], linewidth = 0.58,
        arrow = arrow(length = unit(0.055, "in"), type = "closed")
      ) +
      annotate(
        ggrepel::GeomTextRepel,
        x = 4.48, y = -1.62, label = "lost",
        colour = candidate_transition_colours[["lost"]],
        size = 5.1 / ggplot2::.pt, family = BASE_FAMILY,
        bg.color = "white", bg.r = 0.10,
        max.iter = 0, max.overlaps = Inf,
        box.padding = 0, point.padding = 0, segment.colour = NA
      )
  }
  p
}

rewiring_go <- fread(paths[["rewiring_go"]])
go_similarity <- readRDS(file.path(go_resource_root, "display_wang_similarity.rds"))
semantic_select <- function(
  x,
  n = 5L,
  cutoff = 0.50,
  order_col = "p_adjust",
  overlap_col = "overlap_genes"
) {
  setorderv(x, c(order_col, overlap_col, "go_id"), c(1L, -1L, 1L))
  keep <- integer()
  for (i in seq_len(min(nrow(x), 300L))) {
    if (length(keep)) {
      stopifnot(x$go_id[i] %in% rownames(go_similarity),
                all(x$go_id[keep] %in% colnames(go_similarity)))
      similarity <- go_similarity[x$go_id[i], x$go_id[keep], drop = FALSE]
      if (any(similarity >= cutoff, na.rm = TRUE)) next
    }
    keep <- c(keep, i)
    if (length(keep) >= n) break
  }
  x[keep]
}

hmga1_go <- copy(rewiring_go[
  regulator == "HMGA1" & transition_order == 4L &
    topology_class == "lost" & p_adjust < 0.05
])
hmga1_go <- semantic_select(hmga1_go)
hmga1_go[, term_short := fcase(
  term == "cytoplasmic translation", "Cytoplasmic translation",
  term == "ribonucleoprotein complex biogenesis", "RNP complex biogenesis",
  term == "protein-RNA complex organization", "Protein-RNA organization",
  term == "aerobic electron transport chain", "Aerobic electron transport",
  term == "oxidative phosphorylation", "Oxidative phosphorylation",
  default = term
)]
hmga1_wrapped <- hmga1_go$term_short
hmga1_go[, `:=`(
  gene_ratio = overlap_genes / annotated_module_genes,
  log10_q = -log10(p_adjust),
  term_display = factor(hmga1_wrapped, levels = rev(hmga1_wrapped))
)]

bnc2_gain <- copy(rewiring_go[
  regulator == "BNC2" & transition_order == 3L &
    topology_class == "gained" & p_adjust < 0.05
])
bnc2_loss <- copy(rewiring_go[
  regulator == "BNC2" & transition_order == 4L &
    topology_class == "lost" & p_adjust < 0.05
])
bnc2_common <- merge(
  bnc2_gain,
  bnc2_loss,
  by = "go_id",
  suffixes = c("_gain", "_loss")
)
bnc2_common[, worst_q := pmax(p_adjust_gain, p_adjust_loss)]
bnc2_selected <- semantic_select(
  bnc2_common,
  order_col = "worst_q",
  overlap_col = "overlap_genes_gain"
)
bnc2_go <- rbindlist(list(
  bnc2_selected[, .(
    go_id,
    term = term_gain,
    transition_class = "ProE → Ery\ngained",
    topology_class = "gained",
    overlap_genes = overlap_genes_gain,
    annotated_module_genes = annotated_module_genes_gain
  )],
  bnc2_selected[, .(
    go_id,
    term = term_loss,
    transition_class = "Ery → Norm\nlost",
    topology_class = "lost",
    overlap_genes = overlap_genes_loss,
    annotated_module_genes = annotated_module_genes_loss
  )]
))
bnc2_go[, term_short := fcase(
  term == "immune effector process", "Immune effector process",
  term == "adaptive immune response", "Adaptive immune response",
  term == "leukocyte activation involved in immune response", "Leukocyte activation",
  term == "immune response-activating cell surface receptor signaling pathway",
  "Immune receptor signalling",
  term == "T cell activation", "T-cell activation",
  default = term
)]
bnc2_term_levels <- rev(unique(bnc2_go[topology_class == "gained", term_short]))
bnc2_go[, `:=`(
  transition_class = factor(
    transition_class,
    levels = c("ProE → Ery\ngained", "Ery → Norm\nlost")
  ),
  term_display = factor(term_short, levels = bnc2_term_levels),
  gene_ratio = overlap_genes / annotated_module_genes
)]

candidate_size_scale <- scale_size_continuous(
  name = "Gene ratio",
  range = c(1.8, 4.4),
  limits = c(0, 0.17),
  breaks = c(0.04, 0.08, 0.12),
  labels = percent_format(accuracy = 1),
  guide = guide_legend(
    title.position = "left",
    title.vjust = 0.5,
    nrow = 1,
    byrow = TRUE,
    override.aes = list(fill = "#D8DDE0")
  )
)

p_hmga1_go <- ggplot(hmga1_go, aes(log10_q, term_display)) +
  geom_vline(
    xintercept = -log10(0.05), colour = "#7A858B",
    linetype = "22", linewidth = 0.35
  ) +
  geom_segment(
    aes(x = 0, xend = log10_q, yend = term_display),
    colour = "#D9E0E3", linewidth = 0.35
  ) +
  geom_point(
    aes(size = gene_ratio),
    shape = 21,
    fill = candidate_transition_colours[["lost"]],
    colour = INK_COLOUR, stroke = 0.30
  ) +
  candidate_size_scale +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.13))) +
  labs(x = expression(Ery %->% Norm ~ lost ~ ~ -log[10](BH ~ q)), y = NULL) +
  theme_genome(6.4) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(size = 5.4),
    axis.text.y = element_text(size = 5.2, lineheight = 0.90),
    axis.title.x = element_text(size = 5.6),
    legend.position = "bottom",
    legend.title = element_text(size = 5.4),
    legend.text = element_text(size = 5.1),
    legend.key.size = unit(0.11, "in"),
    legend.margin = margin(0),
    legend.box.margin = margin(0),
    plot.margin = margin(1, 5, 4, 3)
  )

p_bnc2_go <- ggplot(bnc2_go, aes(transition_class, term_display)) +
  geom_point(
    aes(size = gene_ratio, fill = topology_class),
    shape = 21, colour = INK_COLOUR, stroke = 0.30
  ) +
  scale_fill_manual(values = candidate_transition_colours, guide = "none") +
  candidate_size_scale +
  labs(x = NULL, y = NULL) +
  theme_genome(6.4) +
  theme(
    panel.grid.major = element_line(colour = GRID_COLOUR, linewidth = 0.25),
    axis.text.x = element_text(size = 5.2, lineheight = 0.88),
    axis.text.y = element_text(size = 5.2, lineheight = 0.90),
    legend.position = "bottom",
    legend.title = element_text(size = 5.4),
    legend.text = element_text(size = 5.1),
    legend.key.size = unit(0.11, "in"),
    legend.margin = margin(0),
    legend.box.margin = margin(0),
    plot.margin = margin(1, 4, 4, 3)
  )

fig5i <- make_candidate_trajectory("HMGA1") | p_hmga1_go
fig5i <- fig5i +
  patchwork::plot_layout(widths = c(0.95, 1.20)) +
  patchwork::plot_annotation(
    title = "HMGA1",
    theme = theme(
      plot.title = element_text(
        family = BASE_FAMILY, face = "italic", size = 7.4,
        hjust = 0.5, margin = margin(b = 1.5)
      )
    )
  )
fig5j <- make_candidate_trajectory("BNC2") | p_bnc2_go
fig5j <- fig5j +
  patchwork::plot_layout(widths = c(1.03, 0.97)) +
  patchwork::plot_annotation(
    title = "BNC2",
    theme = theme(
      plot.title = element_text(
        family = BASE_FAMILY, face = "italic", size = 7.4,
        hjust = 0.5, margin = margin(b = 1.5)
      )
    )
  )

fig5i <- fig5i & panel_axis_theme(0.8906)
fig5j <- fig5j & panel_axis_theme(0.8567)

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
save_grid_pdf(
  function() {
    grid::grid.newpage()
    grid::grid.draw(fig5d)
  },
  file.path(figure_root, "fig5d.pdf"),
  10,
  6
)
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
save_pdf(fig5i, file.path(figure_root, "fig5i.pdf"), 3.7, 1.50)
save_pdf(fig5j, file.path(figure_root, "fig5j.pdf"), 3.7, 1.50)
message("Saved Fig. 5a-j; Supplementary Fig. S10a-b are rendered by figs10.R.")

source(file.path(repo, "code", "functions", "pdf_assembly.R"), local = TRUE)
assemble_pdf_figure(
  stem = "fig5",
  page = c(183, 211),
  panels = list(
    pdf_panel("a", 10, 200, 36),
    pdf_panel("c", 99, 200, 70),
    pdf_panel("d", 8, 160, 116),
    pdf_panel("e", 127, 156, 47, label_y = 160),
    pdf_panel("f", 8, 114, 74),
    pdf_panel("g", 88, 114, 40, crop = FALSE),
    pdf_panel("h", 131, 114, 40, crop = FALSE),
    pdf_panel("i", 8, 76, 82),
    pdf_panel("j", 93, 76, 81),
    # Draw the network last to keep its lower nodes above panel backgrounds.
    pdf_panel("b", 50, 200, 45)
  ),
  repo = repo,
  figure_root = figure_root
)
