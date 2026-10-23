#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# ---- Parameters ----
network_root <- normalizePath("results/gse162170/network", mustWork = TRUE)
out_dir <- normalizePath("results/gse162170/cell_programs", mustWork = TRUE)
triplet_path <- file.path(network_root, "triplets.tsv")
class_levels <- c("GluN", "IN", "IPC_RG", "OPC_Glia", "Vascular")

target_anchor <- fread(file.path(out_dir, "target_membership.tsv"))
triplets <- fread(triplet_path, showProgress = FALSE)

universe <- unique(triplets$target)
go_rows <- lapply(class_levels, function(anchor) {
  genes <- unique(target_anchor[chromatin_anchor == anchor, target])
  enrich <- enrichGO(
    gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1, readable = FALSE
  )
  x <- as.data.table(as.data.frame(enrich))
  if (!nrow(x)) {
    return(NULL)
  }
  x[, anchor_class := anchor]
  x
})
go_all <- rbindlist(go_rows, use.names = TRUE, fill = TRUE)
go_all[, BH_global := p.adjust(pvalue, method = "BH")]
setnames(go_all, "p.adjust", "BH_within_class")
setorder(go_all, anchor_class, BH_within_class, pvalue, -Count, ID)
go_all[, class_rank := seq_len(.N), by = anchor_class]
go_all[, significant_within_class := BH_within_class < 0.05]
go_all[, display_status := fifelse(significant_within_class, "BH < 0.05", "Not significant")]
go_display <- go_all[,
  {
    if (any(significant_within_class)) {
      .SD[significant_within_class][seq_len(min(3L, sum(significant_within_class)))]
    } else {
      .SD[seq_len(min(2L, .N))]
    }
  },
  by = anchor_class
]
# Significance is a computed result, not a required historical answer.
go_cols <- c(
  "anchor_class", "class_rank", "ID", "Description", "GeneRatio", "BgRatio",
  "pvalue", "BH_within_class", "BH_global", "Count",
  "significant_within_class", "display_status"
)
write_tsv(go_all[, ..go_cols], file.path(out_dir, "go_all_terms.tsv.gz"))
write_tsv(go_display[, ..go_cols], file.path(out_dir, "go_display_terms.tsv"))
