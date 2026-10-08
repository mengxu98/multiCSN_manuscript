#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)
figure_root <- file.path(repo, "figures")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

membership_path <- file.path(
  repo,
  "results",
  "gse162170",
  "cell_programs",
  "module_membership.tsv.gz"
)
summary_path <- file.path(
  repo,
  "results",
  "gse162170",
  "cell_programs",
  "module_summary.tsv"
)
assert_files(c(membership_path, summary_path))

class_order <- c("GluN", "IN", "IPC_RG", "OPC_Glia", "Vascular")
class_labels <- c(
  GluN = "GluN",
  IN = "IN",
  IPC_RG = "IPC/RG",
  OPC_Glia = "OPC/Glia",
  Vascular = "Vascular"
)
class_colours <- c(
  GluN = "#3B6FB6",
  IN = "#8A5BA7",
  IPC_RG = "#D98A2B",
  OPC_Glia = "#4AA684",
  Vascular = "#C4474D"
)
class_annotation_colours <- setNames(
  unname(class_colours[class_order]),
  unname(class_labels[class_order])
)

membership <- unique(fread(membership_path)[, .(
  module_id,
  module_rank,
  regulator,
  anchor_class,
  sign,
  target
)])
module_summary <- fread(summary_path)
stopifnot(
  uniqueN(membership$module_id) == 25L,
  setequal(unique(membership$anchor_class), class_order),
  uniqueN(membership$target) == 5554L
)

module_meta <- unique(membership[, .(
  module_id,
  module_rank,
  regulator,
  anchor_class,
  sign
)])
module_meta <- merge(
  module_meta,
  module_summary[, .(module_id, source_rank = rank)],
  by = "module_id",
  all.x = TRUE,
  sort = FALSE
)
module_meta[is.na(source_rank), source_rank := module_rank]
module_meta[, anchor_class := factor(anchor_class, levels = class_order)]
setorder(module_meta, anchor_class, source_rank, regulator)

target_sets <- split(membership$target, membership$module_id)
module_ids <- module_meta$module_id
n_modules <- length(module_ids)
target_universe <- 11017L

jaccard <- matrix(
  0,
  nrow = n_modules,
  ncol = n_modules,
  dimnames = list(module_ids, module_ids)
)
pair_rows <- vector("list", choose(n_modules, 2L))
pair_index <- 0L
diag(jaccard) <- 1

for (i in seq_len(n_modules - 1L)) {
  set_i <- target_sets[[module_ids[[i]]]]
  for (j in seq.int(i + 1L, n_modules)) {
    set_j <- target_sets[[module_ids[[j]]]]
    overlap <- length(intersect(set_i, set_j))
    union_size <- length(union(set_i, set_j))
    score <- overlap / union_size
    p_value <- phyper(
      overlap - 1L,
      length(set_i),
      target_universe - length(set_i),
      length(set_j),
      lower.tail = FALSE
    )
    jaccard[i, j] <- jaccard[j, i] <- score
    pair_index <- pair_index + 1L
    pair_rows[[pair_index]] <- data.table(
      module_1 = module_ids[[i]],
      module_2 = module_ids[[j]],
      jaccard = score,
      p_value = p_value
    )
  }
}

pair_stats <- rbindlist(pair_rows)
pair_stats[, q_value := p.adjust(p_value, method = "BH")]
q_value <- matrix(
  NA_real_,
  nrow = n_modules,
  ncol = n_modules,
  dimnames = list(module_ids, module_ids)
)
for (i in seq_len(nrow(pair_stats))) {
  row_index <- match(pair_stats$module_1[[i]], module_ids)
  column_index <- match(pair_stats$module_2[[i]], module_ids)
  q_value[row_index, column_index] <- q_value[column_index, row_index] <-
    pair_stats$q_value[[i]]
}

anchor_lookup <- setNames(as.character(module_meta$anchor_class), module_ids)
pair_stats[, same_anchor := anchor_lookup[module_1] == anchor_lookup[module_2]]
stopifnot(
  nrow(pair_stats) == 300L,
  pair_stats[same_anchor == TRUE, .N] == 50L,
  pair_stats[same_anchor == FALSE, .N] == 250L,
  pair_stats[same_anchor == TRUE & q_value < 0.05, .N] == 50L,
  pair_stats[same_anchor == FALSE & q_value < 0.05, .N] == 75L
)

within_class_order <- lapply(class_order, function(anchor) {
  ids <- module_meta[anchor_class == anchor][order(source_rank)]$module_id
  distance_matrix <- 1 - jaccard[ids, ids]
  distance_matrix[distance_matrix < 0] <- 0
  diag(distance_matrix) <- 0
  ids[hclust(as.dist(distance_matrix), method = "average")$order]
})
ordered_ids <- unlist(within_class_order, use.names = FALSE)
ordered_meta <- module_meta[match(ordered_ids, module_id)]
display_labels <- paste0(
  ordered_meta$regulator,
  " (",
  ifelse(ordered_meta$sign == "positive", "+", "-"),
  ")"
)
anchor_factor <- factor(
  class_labels[as.character(ordered_meta$anchor_class)],
  levels = unname(class_labels[class_order])
)
label_colours <- unname(
  class_colours[as.character(ordered_meta$anchor_class)]
)

jaccard_display <- jaccard[ordered_ids, ordered_ids]
q_for_display <- pmax(
  q_value[ordered_ids, ordered_ids],
  .Machine$double.xmin,
  na.rm = FALSE
)
significance_display <- pmin(-log10(q_for_display), 20)
diag(jaccard_display) <- NA_real_
diag(significance_display) <- NA_real_
rownames(jaccard_display) <- colnames(jaccard_display) <- display_labels
rownames(significance_display) <- colnames(
  significance_display
) <- display_labels

font_family <- "Arial"
label_gp <- gpar(
  fontfamily = font_family,
  fontsize = 5.4,
  fontface = "plain",
  col = label_colours
)
small_gp <- gpar(fontfamily = font_family, fontsize = 5.5, fontface = "plain")
title_gp <- gpar(fontfamily = font_family, fontsize = 7.2, fontface = "bold")

overlap_colours <- c(
  low = "#160B39",
  low_mid = "#5B1A82",
  high_mid = "#B5367A",
  high = "#F6D746"
)

make_top_annotation <- function(show_legend) {
  HeatmapAnnotation(
    `Cell class` = anchor_factor,
    col = list(`Cell class` = class_annotation_colours),
    show_annotation_name = FALSE,
    show_legend = show_legend,
    simple_anno_size = unit(0.083, "in"),
    annotation_legend_param = list(
      title_gp = gpar(
        fontfamily = font_family,
        fontsize = 6.2,
        fontface = "plain"
      ),
      labels_gp = small_gp,
      nrow = 1,
      grid_height = unit(0.083, "in"),
      grid_width = unit(0.134, "in")
    )
  )
}

left_annotation <- rowAnnotation(
  `Cell class` = anchor_factor,
  col = list(`Cell class` = class_annotation_colours),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  simple_anno_size = unit(0.083, "in")
)

common_args <- list(
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = anchor_factor,
  column_split = anchor_factor,
  row_gap = unit(0.0315, "in"),
  column_gap = unit(0.0315, "in"),
  width = unit(2.15, "in"),
  height = unit(2.15, "in"),
  border = "#969DA3",
  rect_gp = gpar(col = NA),
  row_names_gp = label_gp,
  column_names_gp = label_gp,
  column_names_rot = 58,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_title = NULL,
  column_title = NULL,
  na_col = "#E6E9EC"
)

figs8a <- do.call(
  Heatmap,
  c(
    list(
      matrix = jaccard_display,
      name = "Jaccard",
      col = colorRamp2(
        c(0, 0.05, 0.15, 0.30),
        unname(overlap_colours)
      ),
      left_annotation = left_annotation,
      top_annotation = make_top_annotation(TRUE),
      heatmap_legend_param = list(
        title = "Jaccard",
        title_gp = small_gp,
        labels_gp = small_gp,
        direction = "horizontal",
        legend_width = unit(1.063, "in"),
        title_position = "topcenter"
      )
    ),
    common_args
  )
)

figs8b <- do.call(
  Heatmap,
  c(
    list(
      matrix = significance_display,
      name = "-log10(q)",
      col = colorRamp2(
        c(0, 2, 5, 20),
        c("#F7F8FA", "#E3DDEB", "#A07BB6", "#542788")
      ),
      top_annotation = make_top_annotation(FALSE),
      show_row_names = TRUE,
      row_names_side = "right",
      heatmap_legend_param = list(
        title = "-log10(BH q)",
        title_gp = small_gp,
        labels_gp = small_gp,
        direction = "horizontal",
        legend_width = unit(1.063, "in"),
        title_position = "topcenter",
        at = c(0, 5, 10, 15, 20),
        labels = c("0", "5", "10", "15", ">=20")
      )
    ),
    common_args[names(common_args) != "show_row_names"]
  )
)

draw_figs8 <- function() {
  draw(
    figs8a + figs8b,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legends = TRUE,
    gap = unit(0.126, "in"),
    padding = unit(c(0.071, 0.079, 0.220, 0.079), "in")
  )
  upViewport(0)
  grid.text(
    "a  Shared target genes",
    x = unit(0.060, "npc"),
    y = unit(0.975, "npc"),
    just = c("left", "top"),
    gp = title_gp
  )
  grid.text(
    "b  Overlap significance",
    x = unit(0.480, "npc"),
    y = unit(0.975, "npc"),
    just = c("left", "top"),
    gp = title_gp
  )
}

save_grid_pdf(draw_figs8, file.path(figure_root, "figs8.pdf"), 5.5, 3.5)

message(
  "Saved Supplementary Fig. S8 module overlap: same-anchor median Jaccard ",
  format(pair_stats[same_anchor == TRUE, median(jaccard)], digits = 3),
  "; cross-anchor median ",
  format(pair_stats[same_anchor == FALSE, median(jaccard)], digits = 3),
  "."
)
