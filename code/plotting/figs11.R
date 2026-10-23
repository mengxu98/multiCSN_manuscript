#!/usr/bin/env Rscript


repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)
source(file.path(repo, "code", "functions", "locus_common.R"), local = TRUE)

# Explicitly set session theme for self-contained execution.
theme_set(theme_genome(7) + theme(text = element_text(family = BASE_FAMILY)))
figure_root <- file.path(repo, "figures")

obj_znf804a <- read_locus_inputs(repo, "ZNF804A")
figs11a <- make_locus_expression(obj_znf804a)
figs11b <- make_locus_key(obj_znf804a)
save_pdf(figs11a, file.path(figure_root, "figs11a.pdf"), 1.5, 2)
save_pdf(figs11b, file.path(figure_root, "figs11b.pdf"), 6.2, 3.4)

obj_kcnc2 <- read_locus_inputs(repo, "KCNC2")
figs11c <- make_locus_track(obj_kcnc2)
figs11d <- make_locus_expression(obj_kcnc2)
figs11e <- make_locus_key(obj_kcnc2)
save_pdf(figs11c, file.path(figure_root, "figs11c.pdf"), 5.2, 2.8)
save_pdf(figs11d, file.path(figure_root, "figs11d.pdf"), 1.5, 2)
save_pdf(figs11e, file.path(figure_root, "figs11e.pdf"), 6.2, 2.3)

source(file.path(repo, "code", "functions", "pdf_assembly.R"), local = TRUE)
assemble_pdf_figure(
  stem = "figs11",
  page = c(165, 190),
  panels = list(
    pdf_panel("a", 8, 186, 42),
    pdf_panel("b", 55, 186, 102),
    pdf_panel("c", 8, 124, 102),
    pdf_panel("d", 116, 124, 41),
    pdf_panel("e", 8, 63, 149)
  ),
  repo = repo,
  figure_root = figure_root
)
message("Saved merged Supplementary Fig. S11 sample-aware candidate-locus panels.")
