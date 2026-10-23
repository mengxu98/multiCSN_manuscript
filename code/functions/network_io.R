source("code/functions/table_io.R")

write_layered_network <- function(fit, output) {
  tf_region_accounting <- fit$accounting$tf_region
  names(tf_region_accounting)[names(tf_region_accounting) == "region"] <- "target"
  tables <- list(
    tf_gene = fit$tf_gene,
    tf_region = fit$tf_region,
    region_gene = fit$region_gene,
    triplets = fit$triplets,
    mediated_tf_gene = fit$mediated_tf_gene,
    tf_gene_accounting = fit$accounting$tf_gene,
    tf_region_accounting = tf_region_accounting,
    region_gene_accounting = fit$accounting$region_gene,
    endpoint_summary = fit$summary
  )
  for (name in names(tables)) {
    write_tsv(tables[[name]], file.path(output, paste0(name, ".tsv")))
  }
  invisible(NULL)
}

read_fit_checkpoint <- function(path) {
  # gzip is available on the supported macOS/Linux hosts; no R.utils dependency.
  if (grepl("[.]gz$", path)) {
    if (!nzchar(Sys.which("gzip"))) stop("gzip is required for checkpoint reads.")
    return(data.table::fread(cmd = paste("gzip -dc", shQuote(path))))
  }
  data.table::fread(path)
}
