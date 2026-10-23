semantic_select <- function(
  x,
  go_similarity,
  n = 5L,
  cutoff = 0.50,
  order_col = "p_adjust",
  overlap_col = "overlap_genes"
) {
  setorderv(x, c(order_col, overlap_col, "go_id"), c(1L, -1L, 1L))
  keep <- integer()
  for (i in seq_len(min(nrow(x), 300L))) {
    if (length(keep)) {
      similarity <- go_similarity[x$go_id[i], x$go_id[keep], drop = FALSE]
      if (any(similarity >= cutoff, na.rm = TRUE)) next
    }
    keep <- c(keep, i)
    if (length(keep) >= n) break
  }
  x[keep]
}
