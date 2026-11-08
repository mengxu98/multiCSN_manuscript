suppressPackageStartupMessages({library(data.table); library(ggplot2); library(grid)})
source("code/functions/common.R")
source("code/functions/pdf_assembly.R")

repo <- normalizePath(getwd(), mustWork = TRUE)
root <- "source_data/supplementary_tables/beeline_pseudotime_20261105"
observations <- fread(file.path(root, "paired_mode_controls.tsv"))
stopifnot(nrow(observations) == 500L, uniqueN(observations$Case) == 500L,
          all(observations$ShuffleFits == 3L))
observations[, Dataset := factor(DatasetType, levels = c("curated", "synthetic"),
                                 labels = c("Curated", "Synthetic"))]
models <- c("GSD", "HSC", "VSC", "mCAD", "dyn-BF", "dyn-BFC", "dyn-CY", "dyn-LI", "dyn-LL", "dyn-TF")
observations[, Model := factor(DatasetGroup, levels = rev(models))]
colours <- c(Curated = "#4D8877", Synthetic = method_colours[["GENIE3"]])
output <- "figures/beeline_s2_compact_20261106"
source_output <- "source_data/supplementary_tables/beeline_s2_compact_20261106"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dir.create(source_output, recursive = TRUE, showWarnings = FALSE)

paired_panel <- function(metric, legend = FALSE) {
  ggplot(observations, aes(.data[[paste0("Static_", metric)]],
                          .data[[paste0("Pseudotime_", metric)]], colour = Dataset, shape = Dataset)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.35, linetype = 2, colour = MUTED_COLOUR) +
    geom_point(size = 0.95, alpha = 0.65, stroke = 0.2) +
    scale_colour_manual(values = colours) +
    scale_shape_manual(values = c(16, 17)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), labels = scales::label_number(accuracy = 0.01)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = scales::label_number(accuracy = 0.01)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = paste("Static", metric), y = paste("Pseudotime", metric)) +
    theme_genome(8) +
    theme(panel.grid.major = element_line(colour = "#EDF0F2", linewidth = 0.2),
          legend.position = if (legend) "inside" else "none",
          legend.position.inside = c(0.77, 0.13),
          legend.title = element_blank(), legend.text = element_text(size = 7.2),
          legend.key.height = unit(9, "pt"), legend.key.width = unit(9, "pt"),
          legend.background = element_blank(), legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(3, 4, 3, 3))
}

limits <- range(c(observations$PseudotimeMinusStatic_AP, observations$PseudotimeMinusShuffled_AP))
limits <- c(floor(limits[1] * 10) / 10, ceiling(limits[2] * 10) / 10)
difference_panel <- function(field, reference, right = FALSE) {
  means <- observations[, .(Difference = mean(get(field))), by = .(Model, Dataset)]
  ggplot(observations, aes(.data[[field]], Model, fill = Dataset)) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = MUTED_COLOUR) +
    geom_boxplot(width = 0.48, outlier.shape = NA, alpha = 0.3, linewidth = 0.3, colour = INK_COLOUR) +
    geom_point(aes(colour = Dataset), size = 0.55, alpha = 0.38,
               position = position_jitter(width = 0, height = 0.12, seed = 18)) +
    geom_point(data = means, aes(Difference, Model), inherit.aes = FALSE,
               shape = 23, fill = MULTICSN_COLOUR, colour = "white", stroke = 0.25, size = 1.8) +
    facet_grid(Dataset ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = colours, guide = "none") +
    scale_colour_manual(values = colours, guide = "none") +
    scale_x_continuous(limits = limits, breaks = seq(-0.2, 1, 0.2),
                       labels = scales::label_number(accuracy = 0.1)) +
    scale_y_discrete(expand = expansion(add = 0.55)) +
    labs(x = paste0("AP difference (pseudotime - ", reference, ")"), y = NULL) +
    theme_genome(8) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "#EDF0F2", linewidth = 0.2),
          panel.spacing.y = unit(2, "mm"),
          strip.text.y = element_text(size = 7.2, colour = if (right) INK_COLOUR else "transparent"),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 7.4, hjust = 0,
                                    colour = if (right) "transparent" else INK_COLOUR),
          axis.title.x = element_text(size = 7.2, margin = margin(t = 5)),
          plot.margin = margin(3, 4, 3, 3), legend.position = "none")
}

panels <- list(a = paired_panel("AP", TRUE), b = paired_panel("AUPRC"),
               c = difference_panel("PseudotimeMinusStatic_AP", "static"),
               d = difference_panel("PseudotimeMinusShuffled_AP", "shuffled", TRUE))
built <- lapply(panels, ggplot_build)
stopifnot(identical(built$a$data[[2]]$x, observations$Static_AP),
          identical(built$a$data[[2]]$y, observations$Pseudotime_AP),
          identical(built$b$data[[2]]$x, observations$Static_AUPRC),
          identical(built$b$data[[2]]$y, observations$Pseudotime_AUPRC),
          identical(sort(built$c$data[[3]]$x), sort(observations$PseudotimeMinusStatic_AP)),
          identical(sort(built$d$data[[3]]$x), sort(observations$PseudotimeMinusShuffled_AP)))
current_boxes <- rbindlist(list(Static = as.data.table(built$c$data[[2]]),
                                Shuffled = as.data.table(built$d$data[[2]])), idcol = "Comparison")
previous_boxes <- fread("source_data/supplementary_tables/beeline_s2_layout_20261105/figs2_box_summaries.tsv")
statistics <- c("Comparison", "group", "xmin", "xlower", "xmiddle", "xupper", "xmax")
setorderv(current_boxes, c("Comparison", "group"))
setorderv(previous_boxes, c("Comparison", "group"))
stopifnot(isTRUE(all.equal(current_boxes[, ..statistics], previous_boxes[, ..statistics],
                          check.attributes = FALSE)))
for (panel_id in names(panels)) {
  height <- if (panel_id %in% c("a", "b")) 59 else 53
  save_pdf(panels[[panel_id]], file.path(output, paste0("figs2", panel_id, ".pdf")), 66 / 25.4, height / 25.4)
  ggsave(file.path(output, paste0("figs2", panel_id, ".svg")), panels[[panel_id]],
         width = 66, height = height, units = "mm", device = svglite::svglite)
  ggsave(file.path(output, paste0("figs2", panel_id, ".png")), panels[[panel_id]],
         width = 66, height = height, units = "mm", dpi = 300, bg = "white")
}
assemble_pdf_figure("figs2", page = c(140, 120), repo = repo, figure_root = output,
                    panels = list(pdf_panel("a", 4, 117, 66, label_x = 3, crop = FALSE),
                                  pdf_panel("b", 73, 117, 66, label_x = 72, crop = FALSE),
                                  pdf_panel("c", 4, 56, 66, label_x = 3, crop = FALSE),
                                  pdf_panel("d", 73, 56, 66, label_x = 72, crop = FALSE)))
stopifnot(system2("pdftoppm", c("-scale-to", "1800", "-singlefile", "-png",
                               file.path(output, "figs2.pdf"), file.path(output, "figs2"))) == 0L,
          system2("pdftocairo", c("-svg", file.path(output, "figs2.pdf"),
                                 file.path(output, "figs2.svg"))) == 0L)
fwrite(current_boxes, file.path(source_output, "figs2_box_summaries.tsv"), sep = "\t")
stopifnot(file.copy(file.path(root, "paired_mode_controls.tsv"),
                    file.path(source_output, "figs2_source.tsv"), overwrite = TRUE))
