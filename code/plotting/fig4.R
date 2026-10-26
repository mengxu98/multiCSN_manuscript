#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(scop)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(scales)
})
data.table::setDTthreads(1L)

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)

# Explicitly set session theme for self-contained execution.
theme_set(theme_genome(7) + theme(text = element_text(family = BASE_FAMILY)))
source(file.path(repo, "code", "functions", "locus_common.R"), local = TRUE)
figure_root <- file.path(repo, "figures")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

source_root <- file.path(
  repo,
  "results", "gse162170", "network"
)
paths <- file.path(
  source_root,
  c(
    "network_topology.tsv",
    "Fig4b_target_funnel.tsv",
    "terminal_states.tsv",
    "chain_summary.tsv"
  )
)

topology <- fread(paths[[1L]])
funnel <- fread(paths[[2L]])
outcomes <- fread(paths[[3L]])
chain_summary <- fread(paths[[4L]])

relation_colours <- c(
  tf_gene = "#3F83B8",
  tf_region = "#8064A2",
  region_gene = "#2A9D8F"
)
topology[, relation_type := factor(relation_type, levels = names(relation_colours))]
topology[,
  relation_set := factor(
    relation_set,
    levels = c("Candidate pairs", "Selected edges")
  )
]
fig4a <- ggplot(
  topology,
  aes(relation_set, relations, group = relation_type, colour = relation_type)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  facet_wrap(
    ~relation_type,
    labeller = as_labeller(c(
      tf_gene = "TF–gene",
      tf_region = "TF–region",
      region_gene = "Region–gene"
    )),
    nrow = 1
  ) +
  scale_colour_manual(values = relation_colours, guide = "none") +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = NULL, y = "Relations (log10)") +
  theme_genome(8) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

funnel[, stage := factor(stage, levels = rev(stage))]
stage_cols <- c(
  "Annotated targets" = "#90A4AE",
  "Genomic domain" = "#78909C",
  "Observed in ≥3 cells" = "#607D8B",
  "Variable candidates" = "#546E7A",
  "Fit eligible" = "#455A64",
  "Selected support" = MULTICSN_COLOUR
)
fig4b <- ggplot(funnel, aes(fraction_of_annotated_targets, stage, fill = stage)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  geom_text(
    aes(
      label = sprintf(
        "%s  (%s)",
        percent(fraction_of_annotated_targets, accuracy = 0.1),
        comma(targets)
      )
    ),
    hjust = 1.02,
    size = 6.2 / ggplot2::.pt,
    family = BASE_FAMILY,
    colour = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(values = stage_cols, guide = "none") +
  scale_x_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, 0.25),
    labels = percent,
    expand = c(0, 0)
  ) +
  labs(x = "Annotated targets", y = NULL) +
  theme_genome(7) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(2, 4, 2, 2)
  )

outcomes[, barrier_label := factor(barrier_label, levels = rev(barrier_label))]
fig4d <- ggplot(outcomes, aes(fraction, barrier_label, fill = barrier_label)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
  geom_text(
    aes(
      label = sprintf(
        "%s  (%s)",
        percent(fraction, accuracy = 0.1),
        comma(targets)
      )
    ),
    hjust = -0.05,
    size = 6.5 / ggplot2::.pt,
    family = BASE_FAMILY,
    colour = INK_COLOUR,
    fontface = "plain"
  ) +
  scale_x_continuous(limits = c(0, 0.70), labels = percent, expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "No eligible domain peak" = "#A8B0B7",
      "Constant RNA response" = "#8064A2",
      "Empty sparse support" = "#E27A24",
      "Selected support" = MULTICSN_COLOUR
    ),
    guide = "none"
  ) +
  labs(x = "Annotated targets", y = NULL) +
  theme_genome(7) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

region_display <- function(x) {
  parts <- tstrsplit(sub("^chr", "", x), "-", fixed = TRUE)
  paste0(
    "chr",
    parts[[1L]],
    ":",
    comma(as.numeric(parts[[2L]])),
    "–",
    comma(as.numeric(parts[[3L]]))
  )
}
znf804a_obj <- read_locus_inputs(repo, "ZNF804A")
display_paths <- copy(znf804a_obj$paths)
setorder(display_paths, start, end, regulator)
display_paths[,
  direction := ifelse(chain_direction >= 0, "Positive", "Negative")
]
display_paths[, backbone := ifelse(backbone_selected, "Present", "Absent")]
region_levels <- unique(display_paths$region)
display_paths[, region_factor := factor(region, levels = region_levels)]
display_paths[,
  regulator_factor := factor(regulator, levels = rev(unique(regulator)))
]
display_paths[,
  direction := factor(direction, levels = c("Negative", "Positive"))
]
display_paths[, backbone := factor(backbone, levels = c("Absent", "Present"))]
fig4j <- ggplot(display_paths, aes(regulator_factor, region_factor)) +
  geom_point(
    aes(fill = direction, shape = backbone),
    size = 1.5,
    colour = INK_COLOUR,
    stroke = 0.45
  ) +
  scale_fill_manual(
    values = c(Negative = "#3B86B6", Positive = MULTICSN_COLOUR)
  ) +
  scale_shape_manual(values = c(Absent = 21, Present = 24)) +
  guides(
    fill = guide_legend(
      order = 1,
      nrow = 1,
      title.position = "left",
      override.aes = list(shape = 21)
    ),
    shape = guide_legend(order = 2, nrow = 1, title.position = "left")
  ) +
  scale_y_discrete(
    labels = setNames(region_display(region_levels), region_levels)
  ) +
  labs(
    x = "Regulators",
    y = "Genomic region",
    fill = "Path sign",
    shape = "Direct TF–gene"
  ) +
  theme_genome(7.5) +
  theme(
    axis.text.x = element_text(
      angle = 72,
      face = "italic", hjust = 1,
      vjust = 1,
      size = 5.2
    ),
    axis.text.y = element_text(size = 5.8),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.spacing.y = unit(-0.05, "in"),
    legend.title = element_text(
      size = 5.5,
      margin = margin(0, 3, 0, 0, "pt")
    ),
    legend.text = element_text(size = 5.2),
    legend.margin = margin(0, 0, 0, 0, "pt"),
    legend.box.margin = margin(-6, 0, 0, 0, "pt"),
    legend.key = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "#E8ECF0")
  )

get_measure <- function(name) as.numeric(chain_summary[measure == name, value])
count_data <- data.table(
  category = c("Exact-region chains", "Direct TF–gene overlap"),
  count = c(get_measure("triplets"), get_measure("backbone_selected"))
)
count_data[, category := factor(category, levels = rev(category))]
fig4c <- ggplot(count_data, aes(count, category, fill = category)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = comma(count)),
    hjust = -0.08,
    family = BASE_FAMILY,
    fontface = "plain",
    size = 6.5 / ggplot2::.pt
  ) +
  scale_x_continuous(
    labels = label_number(scale_cut = cut_short_scale()),
    limits = c(0, max(count_data$count) * 1.3),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "Exact-region chains" = MULTICSN_COLOUR,
      "Direct TF–gene overlap" = "#3F83B8"
    ),
    guide = "none"
  ) +
  labs(x = "Relations", y = NULL) +
  theme_genome(7) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

save_pdf(fig4a, file.path(figure_root, "fig4a.pdf"), 2.8, 1.5)
save_pdf(fig4b, file.path(figure_root, "fig4b.pdf"), 2.5, 1.2)
save_pdf(fig4c, file.path(figure_root, "fig4c.pdf"), 2.8, 0.6)
save_pdf(fig4d, file.path(figure_root, "fig4d.pdf"), 3, 0.8)
save_pdf(fig4j, file.path(figure_root, "fig4j.pdf"), 3.25, 1.9)

source_root <- file.path(repo, "results", "gse162170", "cell_programs")
paths <- file.path(
  source_root,
  c(
    "GSE4a_cell_coordinates_and_labels.tsv",
    "aucell_by_cell.tsv.gz",
    "go_all_terms.tsv.gz"
  )
)

class_order <- c("GluN", "IN", "IPC_RG", "OPC_Glia", "Vascular")
class_labels <- c("GluN", "IN", "IPC/RG", "OPC/Glia", "Vascular")
class_colours <- c(
  GluN = "#3B6FB6",
  IN = "#8A5BA7",
  IPC_RG = "#D98A2B",
  OPC_Glia = "#4AA684",
  Vascular = "#C4474D"
)

coords <- fread(paths[[1L]])
n_cells <- nrow(coords)
dummy <- sparseMatrix(
  i = rep(1L, n_cells),
  j = seq_len(n_cells),
  x = 1,
  dims = c(1L, n_cells),
  dimnames = list("dummy", coords$cell)
)
meta <- as.data.frame(
  coords[, .(broad_class, author_cluster_name, Sample.ID, Sample.Age)]
)
rownames(meta) <- coords$cell
srt <- CreateSeuratObject(dummy, meta.data = meta)
embedding <- as.matrix(coords[, .(dimension_1, dimension_2)])
rownames(embedding) <- coords$cell
colnames(embedding) <- c("UMAP_1", "UMAP_2")
srt[["umap"]] <- CreateDimReducObject(
  embeddings = embedding,
  key = "UMAP_",
  assay = DefaultAssay(srt)
)
aucell <- fread(paths[[2L]])
aucell_wide <- data.table::dcast(
  aucell,
  cell ~ target_program,
  value.var = "AUC"
)
setDT(aucell_wide)
auc_matrix <- t(as.matrix(aucell_wide[, class_order, with = FALSE]))
colnames(auc_matrix) <- aucell_wide$cell
rownames(auc_matrix) <- class_order
feature_map <- setNames(class_order, class_labels)
threshold_path <- file.path(
  repo, "results/gse162170/cell_programs/aucell_thresholds.tsv"
)
threshold_table <- fread(threshold_path)
auc_thresholds <- setNames(threshold_table$threshold, threshold_table$target_program)
fig4e <- scop::CellScoringPlot(
  srt,
  method = "AUCell",
  scores = auc_matrix,
  features = feature_map,
  thresholds = auc_thresholds,
  seed = 42,
  group.by = "broad_class",
  reduction = "umap",
  raster = TRUE,
  group.palcolor = class_colours,
  score.palette = "Spectral",
  highlight.mark.type = "ellipse",
  ellipse.level = 0.80
)
save_pdf(fig4e, file.path(figure_root, "fig4e.pdf"), 9, 3.9)

module_path <- file.path(source_root, "module_raw_and_scaled_activity.tsv")
module_activity <- fread(module_path)
module_activity[, sign_symbol := fifelse(sign == "positive", "+", "-")]
module_activity[,
  module_short := sprintf(
    "%s%s [%d/%d]",
    regulator,
    sign_symbol,
    as.integer(targets),
    as.integer(regions)
  )
]
module_order <- unique(module_activity[
  order(match(anchor_class, class_order), module_rank),
  .(module_id, module_short, anchor_class)
])

make_module_matrix <- function(modality_value) {
  current <- copy(module_activity[modality == modality_value])
  current[, module_id := factor(module_id, levels = module_order$module_id)]
  current[, broad_class := factor(broad_class, levels = class_order)]
  wide <- data.table::dcast(
    current,
    broad_class ~ module_id,
    value.var = "scaled_mean_within_module_modality"
  )
  mat <- as.matrix(wide[
    match(class_order, broad_class),
    module_order$module_id,
    with = FALSE
  ])
  rownames(mat) <- class_labels
  colnames(mat) <- module_order$module_short
  mat
}

module_matrix <- rbind(
  make_module_matrix("RNA_target_AUCell"),
  make_module_matrix("ATAC_selected_region_mean")
)
module_row_labels <- rep(class_labels, 2L)
module_row_colours <- rep(unname(class_colours[class_order]), 2L)
module_anchor <- factor(
  module_order$anchor_class,
  levels = class_order,
  labels = class_labels
)
module_column_colours <- unname(class_colours[module_order$anchor_class])
module_colour_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#4F8FC4", "#FFF7BC", "#D7301F")
)
fig4f <- ComplexHeatmap::Heatmap(
  module_matrix,
  name = "Module\nz score",
  col = module_colour_fun,
  width = ncol(module_matrix) * unit(0.100, "in"),
  height = nrow(module_matrix) * unit(0.100, "in"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = factor(
    rep(
      c("RNA target AUCell", "Selected-region ATAC"),
      each = length(class_order)
    ),
    levels = c("RNA target AUCell", "Selected-region ATAC")
  ),
  row_title = c("RNA\nAUCell", "Region\nATAC"),
  column_split = module_anchor,
  row_gap = unit(0.075, "in"),
  column_gap = unit(0.0433, "in"),
  row_title_side = "left",
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 7.3, fontface = "plain"),
  column_title_gp = gpar(fontsize = 7, fontface = "plain"),
  rect_gp = gpar(col = "grey86", lwd = 0.35),
  row_labels = module_row_labels,
  row_names_gp = gpar(fontsize = 7, col = module_row_colours),
  column_labels = italic_gene_labels(colnames(module_matrix)),
  column_names_gp = gpar(fontsize = 5.8, col = module_column_colours),
  column_names_rot = 62,
  column_names_side = "bottom",
  heatmap_legend_param = list(
    at = c(-2, 0, 2),
    title_gp = gpar(fontsize = 7, fontface = "plain"),
    labels_gp = gpar(fontsize = 7),
    legend_height = unit(0.9, "in")
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    row_class <- class_order[((i - 1L) %% length(class_order)) + 1L]
    if (identical(row_class, module_order$anchor_class[[j]])) {
      grid.circle(
        x,
        y,
        r = unit(0.037, "in"),
        gp = gpar(fill = "white", col = "black", lwd = 0.55)
      )
    }
  }
)
save_grid_pdf(
  function() {
    grid.newpage()
    ComplexHeatmap::draw(fig4f, newpage = FALSE, heatmap_legend_side = "right")
  },
  file.path(figure_root, "fig4f.pdf"),
  4.0,
  2.3
)

go <- fread(paths[[3L]])
go_term_labels <- c(
  "GO:0050808" = "Synapse organization",
  "GO:0048812" = "Neuron projection morphogenesis",
  "GO:0120039" = "Cell projection morphogenesis",
  "GO:0050804" = "Chemical synaptic transmission",
  "GO:0099177" = "Trans-synaptic signalling",
  "GO:0098815" = "Excitatory postsynaptic potential",
  "GO:0051276" = "Chromosome organization",
  "GO:0061351" = "Neural precursor proliferation",
  "GO:0051301" = "Cell division",
  "GO:0001523" = "Retinoid metabolism",
  "GO:0006721" = "Terpenoid metabolism",
  "GO:0002443" = "Leukocyte-mediated immunity",
  "GO:0048514" = "Blood vessel morphogenesis",
  "GO:0001525" = "Angiogenesis",
  "GO:0001568" = "Blood vessel development"
)
matched_term <- go$ID %in% names(go_term_labels)
go[matched_term, Description := unname(go_term_labels[ID])]
go_plot <- copy(go[ID %in% names(go_term_labels)])
ratio_parts <- tstrsplit(go_plot$GeneRatio, "/", fixed = TRUE)
go_plot[,
  gene_ratio := as.numeric(ratio_parts[[1L]]) / as.numeric(ratio_parts[[2L]])
]
go_plot[, `:=`(
  cell_type = factor(
    anchor_class,
    levels = class_order,
    labels = class_labels
  ),
  term = factor(
    unname(go_term_labels[ID]),
    levels = unname(go_term_labels)
  )
)]
fig4g <- ggplot(
  go_plot,
  aes(cell_type, term, size = gene_ratio, fill = BH_within_class)
) +
  geom_point(shape = 21, colour = INK_COLOUR, stroke = 0.32) +
  geom_point(data = go_plot[BH_within_class >= 0.05],
             aes(cell_type, term), inherit.aes = FALSE,
             shape = 4, size = 0.9, colour = INK_COLOUR, stroke = 0.3) +
  scale_size_area(
    name = "Gene ratio",
    max_size = 3.6,
    breaks = c(0.05, 0.10)
  ) +
  scale_fill_gradientn(
    name = "BH-adjusted P",
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    colours = c(MULTICSN_COLOUR, "#F6E8A6", "#A7BBC8"),
    guide = guide_colourbar(
      display = "rectangles",
      frame.colour = INK_COLOUR,
      ticks.colour = INK_COLOUR,
      title = "BH-adjusted P",
      title.position = "top",
      title.hjust = 0,
      direction = "vertical",
      order = 2,
      barwidth = unit(0.075, "in"),
      barheight = unit(0.45, "in")
    )
  ) +
  labs(x = NULL, y = NULL) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "top",
      nrow = 2,
      override.aes = list(fill = "grey55")
    )
  ) +
  theme_genome(7) +
  theme(
    text = element_text(family = BASE_FAMILY),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      vjust = 1,
      size = 6.1,
      colour = INK_COLOUR
    ),
    axis.text.y = element_text(size = 5.8, colour = INK_COLOUR),
    panel.border = element_rect(
      fill = NA,
      colour = "#7A858B",
      linewidth = 0.35
    ),
    panel.grid.major = element_line(colour = "#D9DEE1", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.spacing.y = unit(0.02, "in"),
    legend.key.size = unit(0.12, "in"),
    legend.title = element_text(size = 5.5),
    legend.text = element_text(size = 5.2),
    legend.margin = margin(0, 0, 0, 0, "pt"),
    plot.margin = margin(1, 5, 1, 2, "pt")
  )
save_pdf(fig4g, file.path(figure_root, "fig4g.pdf"), 3, 2)

gsea_root <- file.path(repo, "results", "gse162170", "key_gene_gsea")
gsea_paths <- file.path(
  gsea_root,
  c(
    "gsea.rds",
    "gene_rank_scores.tsv.gz"
  )
)
gsea_result <- readRDS(gsea_paths[[1L]])
gene_ranks <- fread(gsea_paths[[2L]])
# Export the exact ranking and running score used by GSEAPlot.
gsea_curve <- as.data.table(getFromNamespace("gsInfo", "scop")(
  gsea_result$results[["IPC_RG-GO_BP"]], "GO:0006259"
))
gsea_curve <- merge(
  gsea_curve[, .(
    anchor_class = "IPC_RG", pathway = "GO:0006259", gene,
    rank = x, run_score = geneList, running_score = runningScore,
    pathway_member = position == 1L,
    leading_edge = gene %chin% strsplit(CoreGene[[1L]], "/", fixed = TRUE)[[1L]]
  )],
  gene_ranks[anchor_class == "IPC_RG", .(gene, rank_score)],
  by = "gene", sort = FALSE
)
setorder(gsea_curve, rank)
fwrite(gsea_curve, file.path(gsea_root, "gsea_rank_curve.tsv.gz"),
  sep = "\t", quote = FALSE, na = "NA"
)
fig4h <- scop::GSEAPlot(
  res = gsea_result,
  db = "GO_BP",
  plot_type = "line",
  group_use = "IPC_RG",
  id_use = list(IPC_RG = "GO:0006259"),
  padjustCutoff = 0.05,
  direction = "pos",
  n_coregene = 0,
  features_label = c("TOP2A", "CENPF", "PCNA"),
  label.size = 2.35,
  line_width = 1.25,
  line_color = unname(class_colours[["IPC_RG"]]),
  legend.position = "right",
  combine = TRUE,
  verbose = FALSE
)
fig4h <- italicise_gene_grobs(patchwork::patchworkGrob(fig4h), c("TOP2A", "CENPF", "PCNA"))
save_grid_pdf(function() grid::grid.draw(fig4h), file.path(figure_root, "fig4h.pdf"), 4.2, 3.9)

obj <- znf804a_obj
fig4i <- make_locus_track(obj)
save_pdf(fig4i, file.path(figure_root, "fig4i.pdf"), 3, 3)
message("Rendered Fig. 4a-j directly.")

source(file.path(repo, "code", "functions", "pdf_assembly.R"), local = TRUE)
assemble_pdf_figure(
  stem = "fig4",
  page = c(183, 216),
  panels = list(
    pdf_panel("a", 8, 212, 54),
    pdf_panel("b", 69, 212, 51),
    pdf_panel("c", 124, 212, 55),
    pdf_panel("d", 124, 198.5, 55),
    pdf_panel("e", 8, 182.5, 171),
    pdf_panel("f", 8, 109, 85),
    pdf_panel("g", 110, 109, 58),
    pdf_panel("h", 8, 62, 40),
    pdf_panel("i", 52, 62, 43),
    pdf_panel("j", 99, 62, 80)
  ),
  repo = repo,
  figure_root = figure_root
)
