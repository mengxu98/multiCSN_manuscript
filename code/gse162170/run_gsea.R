#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(data.table)
  library(scop)
})
data.table::setDTthreads(1L)

# ---- Parameters ----
rank_path <- normalizePath("results/gse162170/key_gene_gsea/gene_rank_scores.tsv.gz", mustWork = TRUE)
out_dir <- normalizePath("results/gse162170/key_gene_gsea", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

focus_terms <- data.table(
  anchor_class = "IPC_RG",
  pathway = "GO:0006259",
  expected_description = "DNA metabolic process"
)

ranks <- fread(rank_path, showProgress = FALSE)
required_columns <- c("anchor_class", "gene", "rank_score")
stopifnot(
  all(required_columns %in% names(ranks)),
  all(focus_terms$anchor_class %chin% unique(ranks$anchor_class))
)
ranks <- ranks[anchor_class %chin% focus_terms$anchor_class]
stopifnot(
  ranks[, uniqueN(gene), by = anchor_class][, all(V1 == 11017L)],
  !anyNA(ranks[, ..required_columns]),
  all(is.finite(ranks$rank_score))
)

setorder(ranks, anchor_class, -rank_score, gene)
ranks[, run_score := rank_score + rev(seq_len(.N)) * 1e-12,
  by = anchor_class
]

set.seed(20260825)
gsea_result <- scop::RunGSEA(
  geneID = ranks$gene,
  geneScore = ranks$run_score,
  geneID_groups = factor(
    ranks$anchor_class,
    levels = focus_terms$anchor_class
  ),
  IDtype = "symbol",
  result_IDtype = "symbol",
  species = "Homo_sapiens",
  db = "GO_BP",
  minGSSize = 15,
  maxGSSize = 500,
  verbose = TRUE
)

selected <- merge(
  as.data.table(gsea_result$enrichment),
  focus_terms,
  by.x = c("Groups", "ID"),
  by.y = c("anchor_class", "pathway"),
  all = FALSE
)
stopifnot(
  nrow(selected) == nrow(focus_terms),
  all(selected$Description == selected$expected_description),
  all(selected$NES > 0),
  all(is.finite(selected$p.adjust) & selected$p.adjust < 0.05)
)
selected[, core_gene_n := lengths(strsplit(core_enrichment, "/", fixed = TRUE))]
selected[, display_order := match(Groups, focus_terms$anchor_class)]
setorder(selected, display_order)
selected[, display_order := NULL]

saveRDS(gsea_result, file.path(out_dir, "gsea.rds"))
fwrite(
  selected[, .(
    anchor_class = Groups, pathway = ID, Description, setSize,
    enrichmentScore, NES, pvalue, padj = p.adjust, rank,
    leading_edge, core_gene_n, core_enrichment, Database, Version
  )],
  file.path(out_dir, "gsea_selected_terms.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
message("Saved GSEA result for IPC/RG: ", out_dir)
