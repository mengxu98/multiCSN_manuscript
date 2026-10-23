source("code/functions/table_io.R")
source("code/gse274113/evaluate_chain_response.R")

selected_chains_with_directions <- function(root, method, slug, regulators) {
  load_endpoint <- function(name, keys) {
    path <- file.path(root, method, slug, paste0(name, ".tsv.gz"))
    table <- data.table::fread(path)
    required <- c(keys, "direction_sign")
    if (!all(required %in% names(table)) || anyDuplicated(table[, keys, with = FALSE]) ||
        anyNA(table[, keys, with = FALSE])) {
      stop("Invalid selected endpoint: ", path)
    }
    direction <- table$direction_sign
    if (any(!is.na(direction) & (!is.finite(direction) | !direction %in% c(-1, 0, 1)))) {
      stop("Invalid endpoint directions: ", path)
    }
    table
  }
  chains <- load_endpoint("triplet", c("regulator", "region", "target"))
  chains <- chains[regulator %in% regulators,
                   .(regulator, region, target, chain_direction = direction_sign, rank_magnitude)]
  tf_region <- load_endpoint("tf_region", c("regulator", "region"))[
    , .(regulator, region, tf_region_direction = direction_sign)
  ]
  region_gene <- load_endpoint("region_gene", c("region", "target"))[
    , .(region, target, region_gene_direction = direction_sign)
  ]
  if (nrow(merge(unique(chains[, .(regulator, region)]), tf_region,
                 by = c("regulator", "region"))) != nrow(unique(chains[, .(regulator, region)])) ||
      nrow(merge(unique(chains[, .(region, target)]), region_gene,
                 by = c("region", "target"))) != nrow(unique(chains[, .(region, target)]))) {
    stop("Selected chains contain keys absent from their endpoint projections.")
  }
  chains <- merge(chains, tf_region, by = c("regulator", "region"), all.x = TRUE, sort = FALSE)
  merge(chains, region_gene, by = c("region", "target"), all.x = TRUE, sort = FALSE)
}

joined_endpoint_chains <- function(root, method, slug, regulators) {
  tf_region <- data.table::fread(file.path(root, method, slug, "tf_region.tsv.gz"))[
    regulator %in% regulators, .(regulator, region)
  ]
  region_gene <- data.table::fread(file.path(root, method, slug, "region_gene.tsv.gz"))[
    , .(region, target)
  ]
  chains <- merge(tf_region, region_gene, by = "region", allow.cartesian = TRUE, sort = FALSE)
  chains[, `:=`(chain_direction = NA_real_, tf_region_direction = NA_real_,
                 region_gene_direction = NA_real_)]
  chains
}

compare_chain_responses <- function(clusters) {
  methods <- c("multiCSN", "SCENICplus", "Pando")
  pairs <- utils::combn(methods, 2L, simplify = FALSE)
  comparisons <- lapply(pairs, function(pair) {
    left <- clusters[method == pair[[1]], .(regulator, left = joint_positive_fraction)]
    right <- clusters[method == pair[[2]], .(regulator, right = joint_positive_fraction)]
    paired <- merge(left, right, by = "regulator", sort = TRUE)
    difference <- paired$left - paired$right
    difference <- difference[is.finite(difference)]
    p <- if (length(difference) >= 3L && any(difference != 0)) {
      suppressWarnings(stats::wilcox.test(difference, mu = 0, exact = FALSE)$p.value)
    } else NA_real_
    data.frame(
      endpoint = "selected native triplets", metric = "joint ATAC and RNA significance fraction",
      metric_role = "descriptive evidence consistency; not network recovery",
      method_a = pair[[1]], method_b = pair[[2]],
      cluster_unit = "perturbed_TF_mean_across_timepoints", paired_TF_clusters = length(difference),
      mean_difference_a_minus_b = if (length(difference)) mean(difference) else NA_real_,
      median_difference_a_minus_b = if (length(difference)) median(difference) else NA_real_,
      positive_clusters = sum(difference > 0), negative_clusters = sum(difference < 0),
      tied_clusters = sum(difference == 0), wilcoxon_p = p
    )
  })
  result <- data.table::rbindlist(comparisons)
  result[, holm_p_three_pairs := stats::p.adjust(wilcoxon_p, method = "holm")]
  result
}

evaluate_formal_chain_responses <- function(
  endpoints = "results/gse274113/networks",
  truth = "data/gse274113/heldout_pseudobulk_full",
  output = "results/gse274113/evidence",
  annotations = "results/gse274113/comparison"
) {
  methods <- c("multiCSN", "SCENICplus", "Pando")
  timepoints <- c("day 7", "day 9", "day 11", "day 14")
  slugs <- gsub(" ", "_", timepoints, fixed = TRUE)
  effect_summary <- data.table::fread(file.path(truth, "effect_model_summary.tsv"))
  regulators <- sort(unique(as.character(effect_summary$target)))
  if (length(regulators) != 19L || !setequal(effect_summary$timepoint, timepoints)) {
    stop("Held-out truth failed the 19-TF four-timepoint gate.")
  }
  required <- c(
    unlist(lapply(slugs, function(slug) unlist(lapply(methods, function(method) {
      file.path(endpoints, method, slug, c("triplet.tsv.gz", "tf_region.tsv.gz", "region_gene.tsv.gz"))
    })))),
    unlist(lapply(slugs, function(slug) unlist(lapply(c("atac", "rna"), function(modality) {
      file.path(truth, paste0(modality, "_effects"), paste0(slug, "__", regulators, ".tsv.gz"))
    }))))
  )
  if (any(!file.exists(required))) stop("Missing chain-response inputs: ", paste(required[!file.exists(required)], collapse = ", "))
  destinations <- c(output, annotations)
  destinations <- vapply(destinations, normalizePath, character(1), mustWork = FALSE)
  if (destinations[[1]] == destinations[[2]] ||
      any(startsWith(destinations, paste0(rev(destinations), "/")))) {
    stop("Chain metrics and annotations require separate output directories.")
  }
  for (directory in dirname(destinations)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  stages <- vapply(dirname(destinations), function(parent) tempfile(".chain_response_", tmpdir = parent), character(1))
  for (directory in stages) dir.create(directory)
  on.exit({
    unlink(stages, recursive = TRUE)
  }, add = TRUE)
  metric_rows <- list()
  native_metric_rows <- list()
  annotation_rows <- list()
  cursor <- 0L
  for (i in seq_along(timepoints)) {
    chains <- lapply(methods, function(method) selected_chains_with_directions(
      endpoints, method, slugs[[i]], regulators
    ))
    joined <- lapply(methods, function(method) joined_endpoint_chains(
      endpoints, method, slugs[[i]], regulators
    ))
    for (regulator in regulators) {
      regulator_id <- regulator
      filename <- paste0(slugs[[i]], "__", regulator, ".tsv.gz")
      atac <- data.table::fread(file.path(truth, "atac_effects", filename))
      rna <- data.table::fread(file.path(truth, "rna_effects", filename))
      regions <- unique(unlist(lapply(c(chains, joined), function(x) {
        x[regulator == regulator_id, region]
      })))
      region_map <- if (length(regions)) map_regions_to_effect_truth(regions, atac) else NULL
      for (j in seq_along(methods)) {
        evaluated <- evaluate_selected_chain_response(
          chains[[j]], atac, rna, methods[[j]], regulator, timepoints[[i]],
          region_map = region_map
        )
        joined_evaluation <- evaluate_selected_chain_response(
          joined[[j]], atac, rna, methods[[j]], regulator, timepoints[[i]],
          region_map = region_map
        )
        cursor <- cursor + 1L
        metric_rows[[cursor]] <- joined_evaluation$metrics
        native_metric_rows[[cursor]] <- evaluated$metrics
        annotation_rows[[cursor]] <- evaluated$annotations
      }
    }
    cat("chain_response_timepoint\tPASS\t", timepoints[[i]], "\n", sep = "")
  }
  strata <- data.table::rbindlist(metric_rows, use.names = TRUE, fill = TRUE)
  if (nrow(strata) != 228L || anyDuplicated(strata[, .(method, regulator, timepoint)])) {
    stop("Formal chain-response stratum accounting failed.")
  }
  clusters <- summarize_chain_responses(strata)
  write_tsv(strata, file.path(stages[[1]], "chain_response_stratum_metrics.tsv.gz"))
  write_tsv(clusters, file.path(stages[[1]], "chain_response_TF_cluster_metrics.tsv"))
  native_clusters <- summarize_chain_responses(data.table::rbindlist(native_metric_rows))
  write_tsv(compare_chain_responses(native_clusters), file.path(stages[[1]], "chain_response_method_comparisons.tsv"))
  write_tsv(data.table::rbindlist(annotation_rows, use.names = TRUE, fill = TRUE),
            file.path(stages[[2]], "selected_chain_response_annotations.tsv.gz"))
  for (i in seq_along(stages)) {
    publish_result_files(stages[[i]], destinations[[i]])
  }
  invisible(destinations)
}

if (sys.nframe() == 0L) evaluate_formal_chain_responses()
