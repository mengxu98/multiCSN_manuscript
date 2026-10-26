source("code/functions/table_io.R")

export_factorized_chains <- function(
  input = "runs/gse274113/method_fits/multiCSN",
  output = "results/gse274113/chains"
) {
  timepoints <- c("day 7", "day 9", "day 11", "day 14")
  slugs <- gsub(" ", "_", timepoints, fixed = TRUE)
  if (file.exists(output)) stop("Output already exists: ", output)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(".factorized_chains_", tmpdir = dirname(output))
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  summaries <- list()
  for (i in seq_along(slugs)) {
    slug <- slugs[[i]]
    tf_gene <- data.table::fread(file.path(input, slug, "tf_gene.tsv"))
    tf_region <- data.table::fread(file.path(input, slug, "tf_region.tsv"))
    region_gene <- data.table::fread(file.path(input, slug, "region_gene.tsv"))
    chains <- data.table::as.data.table(multiCSN::factorized_chain_support(
      tf_region, region_gene, tf_gene
    ))
    mediated <- data.table::as.data.table(multiCSN::strongest_mediated_projection(chains))
    destination <- file.path(stage, slug)
    dir.create(destination)
    write_tsv(chains[, c("regulator", "region", "target", "weight",
                        "backbone_selected", "backbone_sign_concordant"), with = FALSE],
              file.path(destination, "triplets.tsv.gz"))
    write_tsv(mediated[, c("regulator", "target", "region", "weight",
                          "strongest_path_ties", "strongest_path_sign_consistent"), with = FALSE],
              file.path(destination, "mediated_tf_gene.tsv.gz"))
    summaries[[i]] <- data.frame(
      timepoint = timepoints[[i]], tf_region_edges = nrow(tf_region),
      region_gene_edges = nrow(region_gene), direct_tf_gene_edges = nrow(tf_gene),
      triplets = nrow(chains), backbone_selected_triplets = sum(chains$backbone_selected),
      backbone_sign_concordant_triplets = sum(chains$backbone_sign_concordant),
      mediated_tf_gene_relations = nrow(mediated)
    )
  }
  write_tsv(data.table::rbindlist(summaries), file.path(stage, "timepoint_summary.tsv"))
  if (!file.rename(stage, output)) stop("Unable to publish factorized chains.")
  invisible(output)
}

if (sys.nframe() == 0L) export_factorized_chains()
