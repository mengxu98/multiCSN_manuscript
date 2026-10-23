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

erythroid_dynamic_panel <- function(repo, heroes, network, trajectory, embedding, state_palette) {
  embedding <- data.table::copy(embedding)
  embedding[, state := factor(state, levels = names(state_palette))]
  expression_subset <- file.path(repo, "data/gse194122/erythroid_expression_subset.rds")
  paths <- c(h5ad = file.path(repo, "data/gse194122/processed_multiome.h5ad"))
  dynamic_features <- unique(c(
    heroes$regulator,
    network[order(-represented_states, -stable_chains, target), target]
  ))
  dynamic_meta <- merge(
    trajectory[is.finite(GEX_pseudotime), .(cell_id, GEX_pseudotime)],
    embedding[, .(cell_id, state)],
    by = "cell_id",
    all.x = TRUE,
    sort = FALSE
  )
  dynamic_matrix <- if (file.exists(expression_subset)) {
    cached <- readRDS(expression_subset)
    as(cached[dynamic_features, dynamic_meta$cell_id, drop = FALSE], "CsparseMatrix")
  } else {
    h5_matrix <- HDF5Array::H5ADMatrix(paths[["h5ad"]])
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
  dynamic_feature_programs[, display_order := dynamic_result$feature_metadata[gene, "index"]]
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
  fig5d <- italicise_gene_grobs(regularise_grob_fonts(attr(dynamic_result$plot, "grobs")$full), heroes$regulator)
  save_grid_pdf(function() { grid::grid.newpage(); grid::grid.draw(fig5d) },
    file.path(repo, "figures/fig5d.pdf"), 10, 6)
  invisible(fig5d)
}

erythroid_program_plot <- function(repo, program_name, title, show_tf = TRUE, base_size = 7) {
  scores <- data.table::fread(file.path(repo, "results/gse194122/cases/program_donor_RNA_scores.tsv"))
  scores <- scores[program == program_name]
  if (!show_tf) scores <- scores[signal == "Target-program RNA"]
  summary <- scores[, .(median = median(RNA_score), lower = quantile(RNA_score, 0.25), upper = quantile(RNA_score, 0.75)), by = .(signal, state_order)]
  colours <- c("TF RNA" = "#D2644D", "Target-program RNA" = "#3676A3")
  ggplot2::ggplot(summary, ggplot2::aes(state_order, median, colour = signal, group = signal)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey82", linewidth = 0.25) +
    ggplot2::geom_line(data = scores[signal == "Target-program RNA"], ggplot2::aes(y = RNA_score, group = donor_id), colour = "#B9CADA", alpha = 0.40, linewidth = 0.28) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = signal), alpha = 0.12, colour = NA, show.legend = FALSE) +
    ggplot2::geom_line(ggplot2::aes(linetype = signal), linewidth = 0.65) +
    ggplot2::geom_point(size = 1.35) +
    ggplot2::scale_colour_manual(values = colours, name = NULL) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_linetype_manual(values = c("TF RNA" = "22", "Target-program RNA" = "solid"), name = NULL) +
    ggplot2::scale_x_continuous(breaks = 1:5, labels = c("HSC", "MK/E", "ProE", "Ery", "Norm")) +
    ggplot2::labs(x = NULL, y = "Standardized RNA score", title = title) +
    theme_genome(base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = base_size, hjust = 0.5),
      axis.text = ggplot2::element_text(size = base_size - 0.7), axis.title.y = ggplot2::element_text(size = base_size - 0.4),
      legend.position = "bottom", legend.text = ggplot2::element_text(size = base_size - 0.7),
      legend.key.width = grid::unit(0.18, "in"), legend.margin = ggplot2::margin(0),
      plot.margin = ggplot2::margin(2, 4, 1, 3), panel.grid.minor = ggplot2::element_blank())
}

erythroid_case_GO <- function(repo) {
  hmga <- data.table::fread(file.path(repo, "results/gse194122/cases/HMGA1_GO_all_tests.tsv.gz"))
  hmga <- hmga[set == "core_8_of_10" & background == "HMGA1_source_targets" & go_id %chin% c("GO:0002181", "GO:0042274")]
  hmga <- hmga[, .(regulator = "HMGA1", go_id, term, overlap, annotated_set = annotated_set_genes,
    annotated_universe = annotated_background_genes, term_genes, p_value, p_adjust,
    test_family = "BH across all eligible terms for the 23-gene recurrent-loss core against fixed HMGA1 source-state targets")]
  sox <- data.table::fread(file.path(repo, "results/gse194122/programs/direct_module_GO_supported.tsv"))
  sox <- sox[module_id == "SOX6 (+)" & go_id %chin% c("GO:0007059", "GO:0140014", "GO:0048821")]
  sox <- sox[, .(regulator = "SOX6", go_id, term, overlap, annotated_set, annotated_universe,
    term_genes, p_value, p_adjust = q_global,
    test_family = "Global BH across all signed direct-target programs and eligible GO terms")]
  selected <- data.table::rbindlist(list(hmga, sox))
  selected[, gene_ratio := overlap / annotated_set]
  selected[, display_order := match(go_id, c("GO:0002181", "GO:0042274", "GO:0007059", "GO:0140014", "GO:0048821"))]
  data.table::setorder(selected, display_order)
  data.table::fwrite(selected, file.path(repo, "source_data/files/fig5ij_GO_display.tsv"), sep = "\t")
  selected
}

erythroid_GO_plot <- function(data, title, fill, breaks) {
  labels <- c(`GO:0002181` = "Cytoplasmic\ntranslation", `GO:0042274` = "Ribosomal small\nsubunit biogenesis",
    `GO:0007059` = "Chromosome\nsegregation", `GO:0140014` = "Mitotic nuclear\ndivision", `GO:0048821` = "Erythrocyte\ndevelopment")
  data[, term_display := factor(labels[go_id], levels = rev(labels[go_id]))]
  ggplot2::ggplot(data, ggplot2::aes(-log10(p_adjust), term_display)) +
    ggplot2::geom_vline(xintercept = -log10(0.05), colour = "#7A858B", linetype = "22", linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = -log10(p_adjust), yend = term_display), colour = "#D9E0E3", linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(size = gene_ratio), shape = 21, fill = fill, colour = INK_COLOUR, stroke = 0.3) +
    ggplot2::scale_size_continuous(name = "Gene ratio", range = c(1.8, 4.0), breaks = breaks,
      labels = scales::percent_format(accuracy = 1), limits = c(0, max(breaks, data$gene_ratio)),
      guide = ggplot2::guide_legend(title.position = "top", nrow = 1, override.aes = list(fill = "#D8DDE0"))) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.18))) +
    ggplot2::labs(x = "-log10(BH q)", y = NULL, title = title) +
    theme_genome(7) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 6.6),
      axis.text.y = ggplot2::element_text(lineheight = 0.95), axis.title.x = ggplot2::element_text(size = 6.7),
      plot.title = ggplot2::element_text(size = 6.6, hjust = 0.5), legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 6.4), legend.text = ggplot2::element_text(size = 6.4),
      legend.key.width = grid::unit(0.12, "in"), legend.key.height = grid::unit(0.09, "in"),
      legend.margin = ggplot2::margin(0), plot.margin = ggplot2::margin(2, 6, 1, 2))
}

erythroid_main_case_panels <- function(repo) {
  go <- erythroid_case_GO(repo)
  titles <- c(HMGA1 = "Translation subset\n(16 genes)", SOX6 = "Red-cell development\nsubset (14 genes)")
  programs <- c(HMGA1 = "HMGA1_recurrent_translation", SOX6 = "SOX6_positive_erythroid")
  plots <- lapply(names(programs), function(tf) {
    trajectory <- erythroid_program_plot(repo, programs[[tf]], titles[[tf]]) +
      ggplot2::labs(y = "RNA z score") +
      ggplot2::scale_colour_manual(values = c("TF RNA" = "#D2644D", "Target-program RNA" = "#3676A3"),
        labels = c("TF RNA" = "TF RNA", "Target-program RNA" = "Target subset"), name = NULL) +
      ggplot2::scale_linetype_manual(values = c("TF RNA" = "22", "Target-program RNA" = "solid"),
        labels = c("TF RNA" = "TF RNA", "Target-program RNA" = "Target subset"), name = NULL) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 6.6), axis.text = ggplot2::element_text(size = 6.6),
        axis.title.y = ggplot2::element_text(size = 6.7), legend.text = ggplot2::element_text(size = 6.4),
        legend.key.width = grid::unit(0.12, "in"), panel.grid.major.x = ggplot2::element_blank())
    enrichment <- erythroid_GO_plot(go[regulator == tf],
      if (tf == "HMGA1") "Recurrent loss core\n(23 genes)" else "Positive direct targets\n(607 annotated genes)",
      if (tf == "HMGA1") "#C65445" else "#3676A3",
      if (tf == "HMGA1") c(0.3, 0.6) else c(0.03, 0.07, 0.11))
    patchwork::wrap_plots(trajectory, enrichment, nrow = 1, widths = c(1, 1.2)) +
      patchwork::plot_annotation(title = tf,
        theme = ggplot2::theme(plot.title = ggplot2::element_text(family = BASE_FAMILY, face = "italic", size = 8, hjust = 0.5)))
  })
  list(i = plots[[1]], j = plots[[2]])
}

assemble_erythroid_main_figure <- function(repo, figure_root) {
  assemble_pdf_figure(stem = "fig5", page = c(183, 211), panels = list(
    pdf_panel("a", 10, 200, 36), pdf_panel("c", 99, 200, 70),
    pdf_panel("d", 8, 160, 116), pdf_panel("e", 127, 156, 47, label_y = 160),
    pdf_panel("f", 8, 114, 74), pdf_panel("g", 88, 114, 40, crop = FALSE),
    pdf_panel("h", 131, 114, 40, crop = FALSE), pdf_panel("i", 8, 76, 82),
    pdf_panel("j", 93, 76, 81), pdf_panel("b", 50, 200, 45)), repo = repo, figure_root = figure_root)
}

