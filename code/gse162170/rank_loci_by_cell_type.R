#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))


# ---- Parameters ----
input <- "results/gse162170/loci/candidate_screen/all_target_locus_screen.tsv"
output <- "results/gse162170/loci/by_class"
# ---------------------------------------------------------------
input <- normalizePath(input, mustWork = TRUE)
if (file.exists(output)) stop("Output already exists: ", output, call. = FALSE)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
stage <- paste0(output, ".staging")
if (file.exists(stage)) stop("Staging path already exists: ", stage, call. = FALSE)
dir.create(stage, recursive = FALSE)

x <- fread(input)
classes <- c("GluN", "IN", "IPC_RG", "OPC_Glia", "Vascular")
eligible <- x[eligible == TRUE & chromatin_anchor %chin% classes]
setorder(
  eligible, chromatin_anchor, -sign_balance, -median_chain_delta_bic,
  -target_margin, -chains, target
)
eligible[, class_rank := seq_len(.N), by = chromatin_anchor]
eligible[, class_label := fifelse(
  chromatin_anchor == "IPC_RG", "IPC/RG",
  fifelse(chromatin_anchor == "OPC_Glia", "OPC/Glia", chromatin_anchor)
)]
summary <- eligible[, .(
  eligible_loci = .N,
  top_target = target[which.min(class_rank)],
  top_target_chains = chains[which.min(class_rank)],
  top_target_regions = regions[which.min(class_rank)]
), by = .(chromatin_anchor, class_label)]
summary[, class_order := match(chromatin_anchor, classes)]
setorder(summary, class_order)
summary[, class_order := NULL]
top10 <- eligible[class_rank <= 10L]
if (!setequal(summary$chromatin_anchor, classes) || nrow(top10) != 50L ||
    any(summary$eligible_loci < 10L) || anyDuplicated(eligible[, .(chromatin_anchor, class_rank)])) {
  stop("Cell-type locus-screen self-check failed.", call. = FALSE)
}
fwrite(eligible, file.path(stage, "eligible_loci_by_class.tsv"), sep = "\t")
fwrite(top10, file.path(stage, "top10_loci_by_class.tsv"), sep = "\t")
fwrite(summary, file.path(stage, "cell_type_summary.tsv"), sep = "\t")
if (!file.rename(stage, output)) stop("Unable to promote cell-type locus screen.")
cat("gse162170_locus_screen_by_class\tPASS\t", output, "\n", sep = "")

