canonical_region <- function(x) {
  x <- as.character(x)
  matched <- regexec("^(.+?)[-:]([0-9]+)-([0-9]+)$", x, perl = TRUE)
  pieces <- regmatches(x, matched)
  vapply(pieces, function(piece) {
    if (length(piece) != 4L) {
      return(NA_character_)
    }
    paste0(piece[[2]], ":", piece[[3]], "-", piece[[4]])
  }, character(1))
}

regions_to_granges <- function(values) {
  values <- canonical_region(values)
  if (anyNA(values)) stop("Unparseable region coordinates.", call. = FALSE)
  matched <- regexec("^(.+?):([0-9]+)-([0-9]+)$", values, perl = TRUE)
  pieces <- regmatches(values, matched)
  GenomicRanges::GRanges(
    vapply(pieces, `[[`, character(1), 2L),
    IRanges::IRanges(
      start = as.integer(vapply(pieces, `[[`, character(1), 3L)),
      end = as.integer(vapply(pieces, `[[`, character(1), 4L))
    )
  )
}

map_regions_to_effect_truth <- function(regions, truth) {
  truth <- data.table::as.data.table(truth)
  required <- c("feature", "logFC", "FDR")
  if (!all(required %in% names(truth)) || !nrow(truth) ||
    anyDuplicated(truth$feature) ||
    any(!is.finite(as.numeric(truth$FDR)))) {
    stop("Invalid ATAC effect truth table.", call. = FALSE)
  }
  query <- data.table::data.table(region = unique(as.character(regions)))
  if (!nrow(query) || anyNA(query$region)) {
    stop("No valid predicted regions supplied.", call. = FALSE)
  }
  query[, canonical := canonical_region(region)]
  truth <- truth[, .(
    feature = as.character(feature),
    logFC = as.numeric(logFC), FDR = as.numeric(FDR)
  )]
  truth[, canonical := canonical_region(feature)]
  if (anyNA(query$canonical) || anyNA(truth$canonical)) {
    stop("Unparseable predicted or truth region.", call. = FALSE)
  }

  query_ranges <- regions_to_granges(query$canonical)
  truth_ranges <- regions_to_granges(truth$canonical)
  hits <- GenomicRanges::findOverlaps(
    query_ranges, truth_ranges,
    ignore.strand = TRUE
  )
  result <- query[, .(region)]
  result[, `:=`(
    atac_truth_feature = NA_character_, atac_logFC = NA_real_,
    atac_FDR = NA_real_, overlap_width = NA_integer_,
    overlap_multiplicity = 0L, exact_coordinate_match = FALSE
  )]
  if (!length(hits)) {
    return(result[])
  }

  candidate <- data.table::data.table(
    query_index = S4Vectors::queryHits(hits),
    truth_index = S4Vectors::subjectHits(hits)
  )
  candidate[, overlap_width := as.integer(IRanges::width(
    GenomicRanges::pintersect(
      query_ranges[query_index], truth_ranges[truth_index],
      ignore.strand = TRUE
    )
  ))]
  candidate[, truth_feature := truth$feature[truth_index]]
  multiplicity <- candidate[, .(
    overlap_multiplicity = .N
  ), by = query_index]
  data.table::setorder(
    candidate, query_index, -overlap_width, truth_feature
  )
  selected <- candidate[, .SD[1L], by = query_index]
  selected <- merge(
    selected, multiplicity,
    by = "query_index", sort = FALSE
  )
  selected[, `:=`(
    region = query$region[query_index],
    atac_truth_feature = truth$feature[truth_index],
    atac_logFC = truth$logFC[truth_index],
    atac_FDR = truth$FDR[truth_index],
    exact_coordinate_match =
      query$canonical[query_index] == truth$canonical[truth_index]
  )]
  selected <- selected[, .(
    region, atac_truth_feature, atac_logFC, atac_FDR,
    overlap_width, overlap_multiplicity, exact_coordinate_match
  )]
  result[selected, on = "region", `:=`(
    atac_truth_feature = i.atac_truth_feature,
    atac_logFC = i.atac_logFC,
    atac_FDR = i.atac_FDR,
    overlap_width = i.overlap_width,
    overlap_multiplicity = i.overlap_multiplicity,
    exact_coordinate_match = i.exact_coordinate_match
  )]
  result[]
}

safe_fraction <- function(numerator, denominator) {
  if (denominator > 0L) numerator / denominator else NA_real_
}

evaluate_selected_chain_response <- function(
  chains, atac_truth, rna_truth, method, regulator, timepoint,
  significance_threshold = 0.05
) {
  chains <- data.table::as.data.table(chains)
  required <- c(
    "regulator", "region", "target", "chain_direction",
    "tf_region_direction", "region_gene_direction"
  )
  if (!all(required %in% names(chains))) {
    stop("Chain table is missing endpoint directions.", call. = FALSE)
  }
  regulator_id <- regulator
  chains <- chains[regulator == regulator_id]
  if (anyDuplicated(chains[, .(regulator, region, target)])) {
    stop("Duplicated selected chain keys.", call. = FALSE)
  }
  rna_truth <- data.table::as.data.table(rna_truth)[, .(
    target = as.character(feature),
    rna_logFC = as.numeric(logFC), rna_FDR = as.numeric(FDR)
  )]
  if (!nrow(rna_truth) || anyDuplicated(rna_truth$target) ||
    any(!is.finite(rna_truth$rna_FDR))) {
    stop("Invalid RNA effect truth table.", call. = FALSE)
  }
  if (!nrow(chains)) {
    metrics <- data.frame(
      method = method, regulator = regulator, timepoint = timepoint,
      native_triplets = 0L, evaluable_triplets = 0L,
      evaluable_regions = 0L, evaluable_targets = 0L,
      exact_coordinate_match_fraction = NA_real_,
      multiple_overlap_fraction = NA_real_,
      atac_positive_triplets = 0L, rna_positive_triplets = 0L,
      joint_positive_triplets = 0L,
      gene_direction_evaluable_triplets = 0L,
      atac_direction_evaluable_triplets = 0L,
      mediated_direction_evaluable_triplets = 0L,
      atac_positive_fraction = NA_real_, rna_positive_fraction = NA_real_,
      joint_positive_fraction = NA_real_,
      joint_gene_direction_accuracy = NA_real_,
      joint_atac_direction_accuracy = NA_real_,
      mediated_direction_consistency = NA_real_, stringsAsFactors = FALSE
    )
    return(list(metrics = metrics, annotations = data.table::data.table()))
  }

  region_map <- map_regions_to_effect_truth(chains$region, atac_truth)
  annotated <- merge(chains, region_map, by = "region", all.x = TRUE, sort = FALSE)
  annotated <- merge(annotated, rna_truth, by = "target", all.x = TRUE, sort = FALSE)
  annotated[, `:=`(
    method = method, perturbation_regulator = regulator,
    timepoint = timepoint,
    atac_evaluable = is.finite(atac_FDR),
    rna_evaluable = is.finite(rna_FDR)
  )]
  annotated[, joint_evaluable := atac_evaluable & rna_evaluable]
  annotated[, `:=`(
    atac_positive = atac_evaluable & atac_FDR <= significance_threshold,
    rna_positive = rna_evaluable & rna_FDR <= significance_threshold
  )]
  annotated[, joint_positive := atac_positive & rna_positive]
  annotated[, `:=`(
    gene_direction_evaluable = joint_positive &
      is.finite(chain_direction) & chain_direction != 0 & rna_logFC != 0,
    atac_direction_evaluable = joint_positive &
      is.finite(tf_region_direction) & tf_region_direction != 0 & atac_logFC != 0,
    mediated_direction_evaluable = joint_positive &
      is.finite(region_gene_direction) & region_gene_direction != 0 &
      atac_logFC != 0 & rna_logFC != 0
  )]
  annotated[, `:=`(
    gene_direction_correct = gene_direction_evaluable &
      sign(rna_logFC) == -sign(chain_direction),
    atac_direction_correct = atac_direction_evaluable &
      sign(atac_logFC) == -sign(tf_region_direction),
    mediated_direction_correct = mediated_direction_evaluable &
      sign(rna_logFC) == sign(atac_logFC) * sign(region_gene_direction)
  )]

  atac_n <- sum(annotated$atac_evaluable)
  rna_n <- sum(annotated$rna_evaluable)
  joint_n <- sum(annotated$joint_evaluable)
  mapped <- annotated$atac_evaluable
  gene_direction_n <- sum(annotated$gene_direction_evaluable)
  atac_direction_n <- sum(annotated$atac_direction_evaluable)
  mediated_direction_n <- sum(annotated$mediated_direction_evaluable)
  metrics <- data.frame(
    method = method, regulator = regulator, timepoint = timepoint,
    native_triplets = nrow(chains), evaluable_triplets = joint_n,
    evaluable_regions = data.table::uniqueN(annotated$region[mapped]),
    evaluable_targets = data.table::uniqueN(
      annotated$target[annotated$rna_evaluable]
    ),
    exact_coordinate_match_fraction = safe_fraction(
      sum(annotated$exact_coordinate_match[mapped]), atac_n
    ),
    multiple_overlap_fraction = safe_fraction(
      sum(annotated$overlap_multiplicity[mapped] > 1L), atac_n
    ),
    atac_positive_triplets = sum(annotated$atac_positive),
    rna_positive_triplets = sum(annotated$rna_positive),
    joint_positive_triplets = sum(annotated$joint_positive),
    gene_direction_evaluable_triplets = gene_direction_n,
    atac_direction_evaluable_triplets = atac_direction_n,
    mediated_direction_evaluable_triplets = mediated_direction_n,
    atac_positive_fraction = safe_fraction(
      sum(annotated$atac_positive), atac_n
    ),
    rna_positive_fraction = safe_fraction(
      sum(annotated$rna_positive), rna_n
    ),
    joint_positive_fraction = safe_fraction(
      sum(annotated$joint_positive), joint_n
    ),
    joint_gene_direction_accuracy = safe_fraction(
      sum(annotated$gene_direction_correct), gene_direction_n
    ),
    joint_atac_direction_accuracy = safe_fraction(
      sum(annotated$atac_direction_correct), atac_direction_n
    ),
    mediated_direction_consistency = safe_fraction(
      sum(annotated$mediated_direction_correct), mediated_direction_n
    ),
    stringsAsFactors = FALSE
  )
  list(metrics = metrics, annotations = annotated[])
}
