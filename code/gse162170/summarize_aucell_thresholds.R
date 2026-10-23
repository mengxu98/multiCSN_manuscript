summarize_aucell_thresholds <- function(root = ".") {
  suppressPackageStartupMessages({
    library(data.table)
    library(AUCell)
  })
  out_dir <- file.path(root, "results/gse162170/cell_programs")
  aucell <- fread(file.path(out_dir, "aucell_by_cell.tsv.gz"))
  programs <- c("GluN", "IN", "IPC_RG", "OPC_Glia", "Vascular")
  wide <- dcast(aucell, cell ~ target_program, value.var = "AUC")
  scores <- t(as.matrix(wide[, programs, with = FALSE]))
  colnames(scores) <- wide$cell
  results <- methods::new(
    "aucellResults",
    SummarizedExperiment::SummarizedExperiment(assays = list(AUC = scores))
  )
  set.seed(42)
  thresholds <- getThresholdSelected(AUCell_exploreThresholds(
    results, nCores = 1, plotHist = FALSE, assignCells = TRUE, verbose = FALSE
  ))
  threshold_table <- data.table(
    target_program = programs,
    threshold = sprintf("%.17g", thresholds[programs])
  )
  aucell[, threshold := unname(thresholds[target_program])]
  summary <- aucell[, .(
    threshold = threshold[[1L]],
    above_threshold = sum(AUC >= threshold),
    cells = .N,
    fraction = mean(AUC >= threshold)
  ), by = .(target_program, broad_class)]
  setorder(summary, target_program, broad_class)
  fwrite(threshold_table, file.path(out_dir, "aucell_thresholds.tsv"),
         sep = "\t", quote = FALSE, na = "NA")
  fwrite(summary, file.path(out_dir, "aucell_threshold_summary.tsv"),
         sep = "\t", quote = FALSE, na = "NA")
  invisible(threshold_table)
}

if (sys.nframe() == 0L) summarize_aucell_thresholds()
