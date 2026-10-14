#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})
data.table::setDTthreads(1L)

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(repo, "code", "functions", "common.R"), local = TRUE)
input_root <- file.path(repo, "results", "gse274113", "replicate_deletion_stability")
output_root <- file.path(repo, "figures")

assemble_supp_fig7 <- function() {
  source(file.path(repo, "code", "functions", "pdf_assembly.R"), local = TRUE)
  assemble_pdf_figure(
    stem = "figs7",
    page = c(183, 156),
    panels = list(
      pdf_panel("a", 11.5, 152, 160),
      pdf_panel("b", 11.5, 108, 160),
      pdf_panel("c", 11.5, 56, 48),
      pdf_panel("d", 63.0, 56, 49.5),
      pdf_panel("e", 116.0, 56, 55.5)
    ),
    repo = repo,
    figure_root = output_root
  )
}

input_paths <- c(
  folds = file.path(input_root, "multilayer_decomposition", "all_fold_metrics.tsv"),
  pairwise = file.path(input_root, "pairwise_stability", "pairwise_fold_metrics.tsv"),
  top100 = file.path(input_root, "summary", "top100_decomposition.tsv"),
  rank_decile = file.path(input_root, "multilayer_decomposition", "rank_decile_decomposition.tsv"),
  summary = file.path(input_root, "summary", "timepoint_stability_summary.tsv")
)
if (any(!file.exists(input_paths))) {
  stop(
    "Missing admitted Source Data: ",
    paste(input_paths[!file.exists(input_paths)], collapse = ", ")
  )
}

folds <- fread(input_paths[["folds"]])
pairwise <- fread(input_paths[["pairwise"]])
top100 <- fread(input_paths[["top100"]])
rank_decile <- fread(input_paths[["rank_decile"]])
summary <- fread(input_paths[["summary"]])

timepoints <- c("day 7", "day 9", "day 11", "day 14")
timepoint_colours <- c(
  `day 7` = "#2C7FB8", `day 9` = "#41B6C4",
  `day 11` = "#7A6BB1", `day 14` = "#D97A2B"
)
expected_folds <- c(`day 7` = 4L, `day 9` = 4L, `day 11` = 3L, `day 14` = 3L)
observed_folds <- folds[, .N, by = timepoint]
top_edge_count <- unique(as.integer(top100$stratum_value))
if (nrow(folds) != 14L || anyDuplicated(folds[, .(timepoint, heldout_replicate)]) ||
  !identical(
    setNames(observed_folds$N, observed_folds$timepoint)[timepoints],
    expected_folds
  ) || nrow(top100) != 4L ||
  nrow(pairwise) != 18L ||
  length(top_edge_count) != 1L ||
  !all(top100$edges == top_edge_count) ||
  !setequal(unique(rank_decile$timepoint), timepoints)) {
  stop("The complete fold-rebuilt figure contract failed.")
}

ink <- "#26343E"
muted <- "#697680"
grid_colour <- "#DCE3E7"
base_family <- "Arial"
theme_paper <- function(base_size = 7.0) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.title = element_text(colour = ink, size = base_size),
      axis.text = element_text(colour = ink, size = base_size - 0.4),
      axis.line = element_line(colour = ink, linewidth = 0.35),
      axis.ticks = element_line(colour = ink, linewidth = 0.35),
      panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.25),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(
        colour = ink, face = "plain", size = base_size
      ),
      legend.title = element_text(size = base_size - 0.2, colour = ink),
      legend.text = element_text(size = base_size - 0.4, colour = ink),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_text(
        face = "plain", size = base_size + 1.8,
        colour = ink
      ),
      plot.margin = margin(2, 3, 2, 2)
    )
}

build_s7_accounting_panel <- function() {
  path <- file.path(
    repo, "results", "gse274113", "factorized_chains",
    "timepoint_summary.tsv"
  )
  assert_files(path)
  chains <- fread(path)
  stopifnot(nrow(chains) == 4L, sum(chains$triplets) == 337377L)
  endpoint_colours <- c(
    `TF–gene` = "#3D7FB1", `TF–region` = "#8064A2",
    `Region–gene` = "#289587", Triplets = MULTICSN_COLOUR
  )
  endpoint_labels <- c(
    "TF–gene" = "TF–gene",
    "TF–region" = "TF–region",
    "Region–gene" = "Region–gene",
    "Triplets" = "Triplets"
  )
  flow <- data.table::melt(
    chains,
    id.vars = "timepoint",
    measure.vars = c(
      "direct_tf_gene_edges", "tf_region_edges", "region_gene_edges",
      "triplets"
    ),
    variable.name = "endpoint", value.name = "relations"
  )
  flow[, `:=`(
    endpoint = factor(
      endpoint,
      levels = c(
        "direct_tf_gene_edges", "tf_region_edges", "region_gene_edges",
        "triplets"
      ),
      labels = c("TF–gene", "TF–region", "Region–gene", "Triplets")
    ),
    timepoint = factor(timepoint, levels = timepoints)
  )]
  ggplot(flow, aes(endpoint, relations, fill = endpoint)) +
    geom_col(width = 0.68) +
    geom_text(
      aes(label = sprintf("%.1fK", relations / 1000)),
      vjust = -0.35, size = 2.1, family = base_family, colour = ink
    ) +
    facet_wrap(~timepoint, nrow = 1) +
    scale_fill_manual(values = endpoint_colours) +
    scale_x_discrete(labels = endpoint_labels) +
    scale_y_continuous(
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      expand = expansion(mult = c(0, 0.16))
    ) +
    labs(x = NULL, y = "Selected relations", fill = NULL) +
    theme_paper(7.0) +
    theme(
      axis.text.x = element_text(size = 6.2, angle = 45, hjust = 1, vjust = 1),
      legend.position = "none",
      strip.text = element_text(face = "plain", size = 7.2, colour = ink),
      plot.margin = margin(2, 4, 2, 2)
    )
}

edge_overlap <- rbind(
  folds[, .(timepoint,
    comparison = "Complete vs deletion fold",
    edge_jaccard
  )],
  pairwise[, .(timepoint, comparison = "Deletion fold pair", edge_jaccard)]
)
edge_overlap[, `:=`(
  timepoint = factor(timepoint, levels = timepoints),
  comparison = factor(comparison, levels = c(
    "Complete vs deletion fold", "Deletion fold pair"
  ))
)]
figs7c <- ggplot(edge_overlap, aes(timepoint, edge_jaccard,
  colour = timepoint,
  shape = comparison,
  group = interaction(timepoint, comparison)
)) +
  geom_boxplot(aes(fill = timepoint),
    width = 0.58, outlier.shape = NA,
    linewidth = 0.35, alpha = 0.20, show.legend = FALSE,
    position = position_dodge(width = 0.64)
  ) +
  geom_point(position = position_jitterdodge(
    jitter.width = 0.12, dodge.width = 0.64
  ), size = 1.25, stroke = 0.45, aes(fill = timepoint)) +
  scale_colour_manual(values = timepoint_colours) +
  scale_fill_manual(values = timepoint_colours) +
  scale_shape_manual(values = c(`Complete vs deletion fold` = 16, `Deletion fold pair` = 21)) +
  scale_y_continuous(
    limits = c(0, 0.52), breaks = seq(0, 0.5, 0.1),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  labs(
    x = NULL, y = "Edge Jaccard", colour = NULL, shape = NULL
  ) +
  guides(colour = "none", fill = "none", shape = guide_legend(nrow = 2)) +
  theme_paper(7.0) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = 1, b = 0),
    legend.key.spacing.y = unit(0.01, "in"),
    legend.key.height = unit(0.06, "in"),
    legend.key.width = unit(0.08, "in")
  )

top_measure_columns <- c(
  "rna_model_fraction", "candidate_retention_given_rna_model",
  "fit_eligibility_given_rna_and_candidate", "selection_given_fit_eligible",
  "unconditional_selection_fraction"
)
top100[, (top_measure_columns) := lapply(.SD, as.numeric),
  .SDcols = top_measure_columns
]
top_long <- data.table::melt(
  top100,
  id.vars = c("timepoint", "stratum_value", "edges"),
  measure.vars = top_measure_columns,
  variable.name = "stage", value.name = "fraction"
)
stage_order <- c(
  "rna_model_fraction", "candidate_retention_given_rna_model",
  "fit_eligibility_given_rna_and_candidate", "selection_given_fit_eligible",
  "unconditional_selection_fraction"
)
stage_labels <- c(
  rna_model_fraction = "RNA model",
  candidate_retention_given_rna_model = "Rebuilt candidate",
  fit_eligibility_given_rna_and_candidate = "Fit eligible",
  selection_given_fit_eligible = "Selected (fit)",
  unconditional_selection_fraction = "Selected (total)"
)
top_long[, `:=`(
  timepoint = factor(timepoint, levels = timepoints),
  stage = factor(stage, levels = stage_order)
)]
figs7d <- ggplot(top_long, aes(stage, fraction,
  group = timepoint,
  colour = timepoint
)) +
  geom_line(linewidth = 0.58, na.rm = TRUE) +
  geom_point(size = 1.6, na.rm = TRUE, shape = 21, stroke = 0.45, aes(fill = timepoint)) +
  scale_colour_manual(values = timepoint_colours) +
  scale_fill_manual(values = timepoint_colours) +
  scale_x_discrete(labels = stage_labels) +
  scale_y_continuous(
    limits = c(0, 1.02), breaks = seq(0, 1, 0.2),
    labels = scales::label_percent(accuracy = 1)
  ) +
  labs(
    x = NULL, y = "Retention fraction", colour = NULL
  ) +
  guides(colour = guide_legend(nrow = 2), fill = "none") +
  theme_paper(7.0) +
  theme(
    axis.text.x = element_text(size = 6.0, angle = 45, hjust = 1, vjust = 1),
    legend.position = "bottom",
    legend.margin = margin(t = 1, b = 0),
    legend.key.width = unit(0.08, "in"),
    legend.key.height = unit(0.06, "in")
  )

decile_plot <- data.table::melt(
  rank_decile,
  id.vars = c("timepoint", "stratum_value", "edges"),
  measure.vars = c(
    "selection_given_fit_eligible", "unconditional_selection_fraction"
  ),
  variable.name = "denominator", value.name = "selection_fraction"
)
decile_plot[, `:=`(
  timepoint = factor(timepoint, levels = timepoints),
  weight_decile = as.integer(stratum_value),
  denominator = factor(
    denominator,
    levels = c(
      "selection_given_fit_eligible",
      "unconditional_selection_fraction"
    ),
    labels = c("Conditional on fit eligibility", "All complete edges")
  )
)]
figs7e <- ggplot(decile_plot, aes(weight_decile, selection_fraction,
  colour = timepoint,
  linetype = denominator,
  group = interaction(timepoint, denominator)
)) +
  geom_line(linewidth = 0.55, na.rm = TRUE) +
  geom_point(size = 1.1, na.rm = TRUE) +
  scale_colour_manual(values = timepoint_colours) +
  scale_linetype_manual(values = c("solid", "22")) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 10)
  ) +
  labs(
    x = "Complete-network weight decile", y = "Selection fraction",
    colour = NULL, linetype = NULL
  ) +
  guides(colour = "none", linetype = guide_legend(nrow = 2)) +
  theme_paper(7.0) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = 1, b = 0),
    legend.key.width = unit(0.11, "in"),
    legend.key.height = unit(0.06, "in")
  )

# ---- panel a: cross-method runtime comparison ---
runtime_path <- file.path(
  repo, "results", "gse274113", "method_comparison", "multiomic_runtime_comparison.tsv"
)
if (!file.exists(runtime_path)) stop("Missing runtime comparison table: ", runtime_path)
runtime <- data.table::fread(runtime_path)
runtime[, seconds := as.numeric(inference_seconds)]
runtime[, minutes := seconds / 60]

runtime[method == "multiCSN" & core_configuration == "8 cores (fork)", clean_label := "multiCSN \u00b7 8 cores"]
runtime[method == "multiCSN" & core_configuration == "single core", clean_label := "multiCSN \u00b7 single core"]
runtime[method == "Pando" & core_configuration == "8 cores (doParallel)", clean_label := "Pando \u00b7 8 cores"]
runtime[method == "Pando" & core_configuration == "single core", clean_label := "Pando \u00b7 single core"]
runtime[method == "LINGER", clean_label := "LINGER \u00b7 8 cores"]
runtime[method == "SCENICplus", clean_label := "SCENIC+ \u00b7 12 threads"]

method_cols <- c(
  "multiCSN \u00b7 8 cores" = "#D70440",
  "multiCSN \u00b7 single core" = "#D70440",
  "Pando \u00b7 8 cores" = "#D95F02",
  "Pando \u00b7 single core" = "#D95F02",
  "LINGER \u00b7 8 cores" = "#0072B2",
  "SCENIC+ \u00b7 12 threads" = "#A35C8F"
)

format_time <- function(s) {
  if (s < 60) return(sprintf("%.0f s", s))
  if (s < 3600) return(sprintf("%.1f min", s / 60))
  return(sprintf("%.1f h", s / 3600))
}
runtime[, time_human := sapply(seconds, format_time)]

base_sec <- runtime[clean_label == "multiCSN \u00b7 single core", seconds]
runtime[, factor_vs_single := seconds / base_sec]
runtime[, rel_tag := fifelse(
  clean_label == "multiCSN \u00b7 single core",
  "1.00\u00d7 (Baseline)",
  sprintf("%.2f\u00d7", factor_vs_single)
)]

order_labels <- c(
  "multiCSN \u00b7 8 cores",
  "multiCSN \u00b7 single core",
  "Pando \u00b7 8 cores",
  "Pando \u00b7 single core",
  "LINGER \u00b7 8 cores",
  "SCENIC+ \u00b7 12 threads"
)
runtime[, label := factor(clean_label, levels = rev(order_labels))]
runtime[, y_num := as.numeric(label)]
y_levels <- levels(runtime$label)
y_cols <- method_cols[y_levels]

x_max <- 5200

figs7a <- ggplot(runtime, aes(y = y_num)) +
  annotate("rect", xmin = 0, xmax = 15, ymin = 0.4, ymax = 7.25, fill = "#EDF7ED", alpha = 0.65) +
  annotate("rect", xmin = 15, xmax = 60, ymin = 0.4, ymax = 7.25, fill = "#F5F8FA", alpha = 0.65) +
  annotate("rect", xmin = 60, xmax = x_max, ymin = 0.4, ymax = 7.25, fill = "#FFF5EB", alpha = 0.65) +
  annotate("rect", xmin = 0.1, xmax = 14.5, ymin = 6.70, ymax = 7.20, fill = "#DCEFDC", alpha = 0.95) +
  annotate("text", x = 3.6, y = 6.95, label = "0 \u2013 15 min", size = 2.3, fontface = "bold", colour = "#1B5E20") +
  annotate("rect", xmin = 15.6, xmax = 58, ymin = 6.70, ymax = 7.20, fill = "#E6ECEF", alpha = 0.95) +
  annotate("text", x = 30, y = 6.95, label = "15 \u2013 60 min", size = 2.3, fontface = "bold", colour = "#37474F") +
  annotate("rect", xmin = 62, xmax = x_max * 0.96, ymin = 6.70, ymax = 7.20, fill = "#FFEAD9", alpha = 0.95) +
  annotate("text", x = 550, y = 6.95, label = "1 \u2013 9 h", size = 2.3, fontface = "bold", colour = "#BF360C") +
  geom_segment(
    aes(x = 0, xend = minutes, y = y_num, yend = y_num, colour = clean_label),
    linewidth = 1.15
  ) +
  geom_point(
    aes(x = minutes, y = y_num, colour = clean_label),
    shape = 21, size = 3.3, stroke = 0.85, fill = "white"
  ) +
  geom_point(
    aes(x = minutes, y = y_num, colour = clean_label),
    shape = 16, size = 2.0
  ) +
  geom_text(
    aes(
      x = minutes * 1.12 + 0.3, y = y_num,
      label = sprintf("%.0f s (%s)  [%s]", seconds, time_human, rel_tag)
    ),
    hjust = 0, size = 2.15, fontface = "bold", colour = ink
  ) +
  scale_y_continuous(
    breaks = 1:6,
    labels = y_levels,
    limits = c(0.4, 7.35),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 1),
    limits = c(0, x_max),
    breaks = c(0, 5, 15, 60, 240, 540),
    labels = c("0 min", "5 min", "15 min", "1 h", "4 h", "9 h"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_colour_manual(values = method_cols, guide = "none") +
  labs(
    x = "Total time",
    y = NULL
  ) +
  theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(
    axis.title = element_text(colour = ink, size = 7.2),
    axis.text.x = element_text(colour = ink, size = 6.8),
    axis.line = element_line(colour = ink, linewidth = 0.35),
    axis.ticks = element_line(colour = ink, linewidth = 0.35),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.6, face = "bold", colour = y_cols),
    panel.grid.major.x = element_line(colour = "#EAEFF2", linewidth = 0.3, linetype = "dotted"),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(2, 6, 2, 2)
  )

save_pdf(figs7a, file.path(output_root, "figs7a.pdf"), 6.30, 1.55)
figs7b <- build_s7_accounting_panel()
save_pdf(figs7b, file.path(output_root, "figs7b.pdf"), 6.30, 1.75)
save_pdf(figs7c, file.path(output_root, "figs7c.pdf"), 1.89, 1.7)
save_pdf(figs7d, file.path(output_root, "figs7d.pdf"), 1.95, 1.7)
save_pdf(figs7e, file.path(output_root, "figs7e.pdf"), 2.19, 1.7)

assemble_supp_fig7()
message("Rendered Supplementary Fig. S7a-e and assembled S7 independently.")
