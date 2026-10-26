#!/usr/bin/env Rscript

source("code/functions/table_io.R")

source("code/functions/ranking_metrics.R")

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})
setDTthreads(4L)

# ---- Parameters ----
# Parameters: NETWORK_ROOT, ENDPOINT_ROOT, TRUTH_ROOT, OUTPUT_ROOT.
arguments <- list(
  network_root = "runs/gse274113/LINGER",
  endpoint_root = "results/gse274113/networks",
  truth_root = "data/gse274113/heldout_pseudobulk_full",
  output_root = "results/gse274113/neighborhood"
)
# ---------------------------------------------------------------

network_root <- normalizePath(arguments$network_root, mustWork = TRUE)
endpoint_root <- normalizePath(arguments$endpoint_root, mustWork = TRUE)
truth_root <- normalizePath(arguments$truth_root, mustWork = TRUE)
output_root <- normalizePath(arguments$output_root, mustWork = FALSE)

timepoints <- c("day 7", "day 9", "day 11", "day 14")
day_slugs <- gsub(" ", "_", timepoints, fixed = TRUE)
names(day_slugs) <- timepoints
targets <- sort(unique(fread(file.path(truth_root, "effect_model_summary.tsv"))$target))

enrichment_row <- function(stats, members, modality, regulator, timepoint, reason) {
  members <- intersect(unique(members), names(stats))
  if (length(members) < 5L || length(members) > length(stats) - 5L) {
    return(data.table(
      pathway = paste("LINGER", modality, "order1", sep = "|"),
      pval = NA_real_, padj = NA_real_, log2err = NA_real_, ES = NA_real_,
      NES = NA_real_, size = length(members), method = "LINGER",
      neighborhood_order = 1L, modality, leading_edge = NA_character_,
      perturbation_regulator = regulator, timepoint, evaluable = FALSE,
      unavailable_reason = if (nzchar(reason)) reason else "pathway_size_out_of_range"
    ))
  }
  result <- suppressWarnings(fgseaMultilevel(
    pathways = list(LINGER = members), stats = stats,
    minSize = 5L, maxSize = length(stats) - 5L, eps = 0,
    scoreType = "pos", nproc = 1L
  ))
  data.table(
    pathway = paste("LINGER", modality, "order1", sep = "|"),
    pval = result$pval[[1]], padj = result$padj[[1]],
    log2err = result$log2err[[1]], ES = result$ES[[1]], NES = result$NES[[1]],
    size = result$size[[1]], method = "LINGER", neighborhood_order = 1L,
    modality, leading_edge = paste(result$leadingEdge[[1]], collapse = ";"),
    perturbation_regulator = regulator, timepoint, evaluable = TRUE,
    unavailable_reason = NA_character_
  )
}

response_row <- function(truth, members, order, regulator, timepoint, reason) {
  members <- intersect(unique(members), truth$feature)
  if (!length(members)) {
    return(data.table(
      set_size = 0L, responsive = 0L, response_rate = NA_real_,
      background_response_rate = NA_real_, response_odds_ratio = NA_real_,
      response_p_value = NA_real_, direction_evaluable = 0L,
      direction_consistent = 0L, direction_consistency = NA_real_,
      method = "LINGER", perturbation_regulator = regulator, timepoint,
      modality = "RNA", neighborhood_order = order, evaluable = FALSE,
      unavailable_reason = reason
    ))
  }
  in_set <- truth[feature %chin% members]
  outside <- truth[!feature %chin% members]
  if (!nrow(outside)) {
    return(data.table(
      set_size = nrow(in_set), responsive = sum(in_set$FDR < 0.05),
      response_rate = mean(in_set$FDR < 0.05), background_response_rate = NA_real_,
      response_odds_ratio = NA_real_, response_p_value = NA_real_,
      direction_evaluable = 0L, direction_consistent = 0L,
      direction_consistency = NA_real_, method = "LINGER",
      perturbation_regulator = regulator, timepoint, modality = "RNA",
      neighborhood_order = order, evaluable = FALSE,
      unavailable_reason = "no_outside_background"
    ))
  }
  in_sig <- sum(in_set$FDR < 0.05)
  out_sig <- sum(outside$FDR < 0.05)
  test <- fisher.test(matrix(
    c(in_sig, nrow(in_set) - in_sig, out_sig, nrow(outside) - out_sig),
    nrow = 2L, byrow = TRUE
  ))
  data.table(
    set_size = nrow(in_set), responsive = in_sig,
    response_rate = in_sig / nrow(in_set),
    background_response_rate = out_sig / nrow(outside),
    response_odds_ratio = unname(test$estimate), response_p_value = test$p.value,
    direction_evaluable = 0L, direction_consistent = 0L,
    direction_consistency = NA_real_, method = "LINGER",
    perturbation_regulator = regulator, timepoint, modality = "RNA",
    neighborhood_order = order, evaluable = is.finite(unname(test$estimate)),
    unavailable_reason = NA_character_
  )
}

enrichment_rows <- list()
response_rows <- list()
for (day_index in seq_along(timepoints)) {
  timepoint <- timepoints[[day_index]]
  slug <- day_slugs[[timepoint]]
  tf_gene <- unique(fread(
    file.path(network_root, slug, "networks", "tf_gene.tsv.gz"),
    select = c("regulator", "target")
  ))
  tf_gene <- tf_gene[target != regulator]
  setkey(tf_gene, regulator)
  regulators <- unique(tf_gene$regulator)

  tf_region <- unique(fread(
    file.path(endpoint_root, "LINGER", slug, "tf_region.tsv.gz"),
    select = c("regulator", "region")
  ))
  tf_region[, feature := sub(":", "-", region, fixed = TRUE)]
  setkey(tf_region, regulator)

  direct_gene <- setNames(vector("list", length(targets)), targets)
  second_gene <- setNames(vector("list", length(targets)), targets)
  direct_region <- setNames(vector("list", length(targets)), targets)
  second_cache <- new.env(parent = emptyenv())
  for (regulator in targets) {
    direct <- tf_gene[regulator, unique(target), on = "regulator"]
    direct <- setdiff(direct, regulator)
    direct_gene[[regulator]] <- direct
    direct_region[[regulator]] <- tf_region[
      regulator, unique(feature), on = "regulator"
    ]
    intermediates <- sort(intersect(direct, regulators))
    cache_key <- paste(intermediates, collapse = ";")
    if (!length(intermediates)) {
      second <- character()
    } else if (exists(cache_key, envir = second_cache, inherits = FALSE)) {
      second <- get(cache_key, envir = second_cache, inherits = FALSE)
    } else {
      second <- tf_gene[regulator %chin% intermediates, unique(target)]
      assign(cache_key, second, envir = second_cache)
    }
    second_gene[[regulator]] <- setdiff(second, c(regulator, direct))
  }
  rm(tf_gene, tf_region, second_cache)
  gc(verbose = FALSE)

  evaluated <- parallel::mclapply(seq_along(targets), function(target_index) {
    regulator <- targets[[target_index]]
    set.seed(274113L + day_index * 100L + target_index)
    rna <- fread(file.path(
      truth_root, "rna_effects", paste0(slug, "__", regulator, ".tsv.gz")
    ))
    atac <- fread(file.path(
      truth_root, "atac_effects", paste0(slug, "__", regulator, ".tsv.gz")
    ))
    rna_reason <- if (length(direct_gene[[regulator]])) "" else
      "perturbed_TF_not_exported_by_LINGER"
    atac_reason <- if (length(direct_region[[regulator]])) "" else
      "perturbed_TF_not_exported_by_LINGER"
    list(
      enrichment = rbindlist(list(
        enrichment_row(
          perturbation_response_stat(rna), direct_gene[[regulator]], "RNA", regulator,
          timepoint, rna_reason
        ),
        enrichment_row(
          perturbation_response_stat(atac), direct_region[[regulator]], "ATAC", regulator,
          timepoint, atac_reason
        )
      )),
      response = rbindlist(list(
        response_row(
          rna, direct_gene[[regulator]], 1L, regulator, timepoint, rna_reason
        ),
        response_row(
          rna, second_gene[[regulator]], 2L, regulator, timepoint,
          if (length(direct_gene[[regulator]]))
            "no_distinct_second_order_targets" else rna_reason
        )
      ))
    )
  }, mc.cores = min(12L, length(targets)), mc.preschedule = FALSE)

  enrichment_rows[[timepoint]] <- rbindlist(lapply(evaluated, `[[`, "enrichment"))
  response_rows[[timepoint]] <- rbindlist(lapply(evaluated, `[[`, "response"))
}

enrichment <- rbindlist(enrichment_rows)
response <- rbindlist(response_rows)
response[, response_p_adjust := NA_real_]
response[is.finite(response_p_value),
  response_p_adjust := p.adjust(response_p_value, method = "BH"),
  by = .(modality, neighborhood_order)
]

dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
stage <- tempfile(".result-files-", tmpdir = dirname(output_root))
dir.create(stage)
on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
fwrite(
  enrichment, file.path(stage, "fig3d_linger_direct_enrichment.tsv.gz"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  response, file.path(stage, "fig3g_linger_neighborhood_response.tsv.gz"),
  sep = "\t", quote = FALSE, na = "NA"
)
publish_result_files(stage, output_root)
