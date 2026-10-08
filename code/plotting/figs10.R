#!/usr/bin/env Rscript


repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)
source(file.path(repo, "code", "functions", "locus_common.R"), local = TRUE)

# Explicitly set session theme for self-contained execution.
theme_set(theme_genome(7) + theme(text = element_text(family = BASE_FAMILY)))
figure_root <- file.path(repo, "figures")

obj_znf804a <- read_locus_inputs(repo, "ZNF804A")
stopifnot(
  nrow(obj_znf804a$paths) == 32L,
  data.table::uniqueN(obj_znf804a$paths$region) == 9L
)
figs10a <- make_locus_expression(obj_znf804a)
figs10b <- make_locus_key(obj_znf804a)
save_pdf(figs10a, file.path(figure_root, "figs10a.pdf"), 1.5, 2)
save_pdf(figs10b, file.path(figure_root, "figs10b.pdf"), 6.2, 3.4)

obj_kcnc2 <- read_locus_inputs(repo, "KCNC2")
stopifnot(
  nrow(obj_kcnc2$paths) == 33L,
  data.table::uniqueN(obj_kcnc2$paths$region) == 9L,
  sum(obj_kcnc2$paths$chain_direction > 0) == 18L,
  sum(obj_kcnc2$paths$chain_direction < 0) == 15L
)
figs10c <- make_locus_track(obj_kcnc2)
figs10d <- make_locus_expression(obj_kcnc2)
figs10e <- make_locus_key(obj_kcnc2)
save_pdf(figs10c, file.path(figure_root, "figs10c.pdf"), 5.2, 2.8)
save_pdf(figs10d, file.path(figure_root, "figs10d.pdf"), 1.5, 2)
save_pdf(figs10e, file.path(figure_root, "figs10e.pdf"), 6.2, 2.3)

source(file.path(repo, "code", "functions", "pdf_assembly.R"), local = TRUE)
assemble_pdf_figure(
  stem = "figs10",
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
message("Saved merged Supplementary Fig. S10 sample-aware candidate-locus panels.")
