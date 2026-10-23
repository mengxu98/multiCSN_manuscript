suppressPackageStartupMessages({
  library(ggplot2)
})

options(Seurat.object.assay.version = "v3")
if (!interactive()) options(device = function(...) grDevices::pdf(file = NULL, ...))

# Register Arial in standard PostScript/PDF device databases to avoid
# grid.Call(C_textBounds) warnings during plot construction
if (!"Arial" %in% names(grDevices::pdfFonts())) {
  grDevices::pdfFonts(Arial = grDevices::pdfFonts()$Helvetica)
}
if (!"Arial" %in% names(grDevices::postscriptFonts())) {
  grDevices::postscriptFonts(Arial = grDevices::postscriptFonts()$Helvetica)
}

MULTICSN_COLOUR <- "#D70440"
INK_COLOUR <- "#24323D"
MUTED_COLOUR <- "#6F7C86"
GRID_COLOUR <- "#DCE3E7"
BASE_FAMILY <- "Arial"

method_colours <- c(
  multiCSN = MULTICSN_COLOUR,
  GENIE3 = "#2177B8",
  GRNBoost2 = "#5E4FA2",
  PIDC = "#2C7FB8",
  PPCOR = "#0AA344",
  LEAP = "#F9BD10",
  SINCERITIES = "#8076A3",
  SCODE = "#CC6F00",
  SCORPION = "#4E9F8A",
  GNIPLR = "#B26D9E",
  SCSGL = "#8A7D55",
  Pando = "#FA8C35",
  LINGER = "#0072B2",
  `SCENIC+` = "#A35C8F",
  SCENICplus = "#A35C8F"
)

theme_genome <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = INK_COLOUR),
      axis.ticks = element_line(linewidth = 0.35, colour = INK_COLOUR),
      axis.text = element_text(colour = INK_COLOUR, size = base_size - 0.3),
      axis.title = element_text(colour = INK_COLOUR, size = base_size),
      panel.grid.major.y = element_line(colour = GRID_COLOUR, linewidth = 0.25),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "plain", colour = INK_COLOUR, size = base_size),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.4),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(3, 4, 3, 4)
    )
}

save_pdf <- function(plot, path, width_in, height_in) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(
    filename = path,
    width = width_in,
    height = height_in,
    family = BASE_FAMILY,
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()
  invisible(path)
}

save_grid_pdf <- function(draw_fun, path, width_in, height_in) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(
    filename = path,
    width = width_in,
    height = height_in,
    family = BASE_FAMILY,
    bg = "white"
  )
  draw_fun()
  grDevices::dev.off()
  invisible(path)
}

italic_gene_labels <- function(labels) {
  as.expression(lapply(as.character(labels), function(label) {
    symbol <- sub("[[:space:]+:].*$", "", label)
    suffix <- substring(label, nchar(symbol) + 1L)
    bquote(italic(.(symbol)) * .(suffix))
  }))
}

italicise_gene_grobs <- function(grob, symbols) {
  if (!inherits(grob, "grob")) return(grob)
  if (inherits(grob, "text") && is.character(grob$label)) {
    matched <- grob$label %in% symbols
    if (any(matched)) {
      faces <- rep_len(if (is.null(grob$gp$font)) 1L else grob$gp$font, length(matched))
      faces[matched] <- ifelse(faces[matched] %in% c(2L, 4L), 4L, 3L)
      grob$gp$font <- faces
    }
  }
  if (inherits(grob, "textrepeltree") && "label" %in% names(grob$data)) {
    matched <- grob$data$label %in% symbols
    grob$data$fontface[matched] <- ifelse(grob$data$fontface[matched] %in% c("bold", "bold.italic"), "bold.italic", "italic")
  }
  if (!is.null(grob$children)) grob$children <- do.call(grid::gList,
    lapply(grob$children, italicise_gene_grobs, symbols = symbols))
  if (!is.null(grob$grobs)) grob$grobs <- lapply(grob$grobs, italicise_gene_grobs, symbols = symbols)
  if (inherits(grob, "recordedGrob") && !is.null(grob$list)) grob$list <- lapply(grob$list,
    function(x) if (inherits(x, "grob")) italicise_gene_grobs(x, symbols) else x)
  grob
}
