#!/usr/bin/env Rscript

source("code/functions/table_io.R")

source("code/functions/ranking_metrics.R")

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
truth_root <- file.path(repo, "data", "gse274113", "heldout_pseudobulk_full")
network_root <- file.path(repo, "results", "gse274113", "networks")
annotation_path <- file.path(repo, "results", "gse274113", "comparison", "selected_chain_response_annotations.tsv.gz")
output <- file.path(repo, "results", "gse274113", "neighborhood")
methods <- c("multiCSN", "Pando", "SCENICplus")
day_slugs <- c("day_7", "day_9", "day_11", "day_14")
day_labels <- c("day 7", "day 9", "day 11", "day 14")
names(day_labels) <- day_slugs


stage <- tempfile(pattern = ".perturbation_neighborhood_", tmpdir = dirname(output))
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
set.seed(274113L)

collapse_tf_gene <- function(path) {
  x <- fread(path, select = c("regulator", "target", "rank_magnitude", "direction_sign", "score_available"))
  x <- x[score_available == TRUE & is.finite(rank_magnitude) & is.finite(direction_sign)]
  if (anyDuplicated(x[, .(regulator, target)])) {
    setorder(x, regulator, target, -rank_magnitude, -direction_sign)
    x <- x[, .SD[1L], by = .(regulator, target)]
  }
  x[, .(regulator, target, rank_magnitude, direction_sign = as.integer(sign(direction_sign)))]
}

build_neighborhood <- function(edges, source_tf) {
  direct <- edges[regulator == source_tf & target != source_tf]
  if (!nrow(direct)) {
    return(data.table())
  }
  direct_out <- direct[, .(
    target,
    order = 1L,
    intermediate = NA_character_,
    path_direction = direction_sign,
    paths = 1L
  )]
  second_edges <- edges[regulator %chin% direct$target]
  if (!nrow(second_edges)) {
    return(direct_out)
  }
  second <- merge(
    direct[, .(intermediate = target, first_direction = direction_sign)],
    second_edges[, .(intermediate = regulator, target, second_direction = direction_sign)],
    by = "intermediate",
    allow.cartesian = TRUE
  )
  second <- second[target != source_tf & !target %chin% direct$target]
  second[, path_direction_each := first_direction * second_direction]
  second_out <- second[, .(
    order = 2L,
    intermediate = paste(sort(unique(intermediate)), collapse = ";"),
    path_direction = if (uniqueN(path_direction_each) == 1L) path_direction_each[1L] else NA_integer_,
    paths = .N
  ), by = target]
  rbindlist(list(direct_out, second_out), use.names = TRUE, fill = TRUE)
}

run_enrichment <- function(stats, pathways, modality, method_names, orders) {
  pathways <- lapply(pathways, intersect, y = names(stats))
  sizes <- lengths(pathways)
  keep <- sizes >= 5L & sizes <= max(5L, length(stats) - 5L)
  if (!any(keep)) {
    return(data.table())
  }
  pathways_use <- pathways[keep]
  result <- suppressWarnings(fgseaMultilevel(
    pathways = pathways_use,
    stats = stats,
    minSize = 5L,
    maxSize = length(stats) - 5L,
    eps = 0,
    scoreType = "pos"
  ))
  result <- as.data.table(result)
  result[, `:=`(
    method = method_names[match(pathway, names(pathways))],
    neighborhood_order = orders[match(pathway, names(pathways))],
    modality = modality,
    leading_edge = vapply(leadingEdge, paste, collapse = ";", character(1))
  )]
  result[, leadingEdge := NULL]
  result[]
}

evaluate_set <- function(truth, members) {
  members <- members[target %chin% truth$feature]
  if (!nrow(members)) {
    return(data.table())
  }
  members <- members[, .SD[1L], by = target]
  joined <- merge(members, truth[, .(target = feature, logFC, FDR)], by = "target")
  in_set_sig <- sum(joined$FDR < 0.05, na.rm = TRUE)
  in_set_non <- nrow(joined) - in_set_sig
  outside <- truth[!feature %chin% joined$target]
  out_sig <- sum(outside$FDR < 0.05, na.rm = TRUE)
  out_non <- nrow(outside) - out_sig
  test <- stats::fisher.test(matrix(c(in_set_sig, in_set_non, out_sig, out_non), nrow = 2L, byrow = TRUE))
  directional <- joined[FDR < 0.05 & is.finite(path_direction)]
  data.table(
    set_size = nrow(joined),
    responsive = in_set_sig,
    response_rate = in_set_sig / nrow(joined),
    background_response_rate = out_sig / nrow(outside),
    response_odds_ratio = unname(test$estimate),
    response_p_value = test$p.value,
    direction_evaluable = nrow(directional),
    direction_consistent = sum(sign(directional$logFC) == -directional$path_direction),
    direction_consistency = if (nrow(directional)) mean(sign(directional$logFC) == -directional$path_direction) else NA_real_
  )
}

effect_summary <- fread(file.path(truth_root, "effect_model_summary.tsv"))
targets <- sort(unique(effect_summary$target))
if (!length(targets) || anyNA(targets) || any(!nzchar(targets))) stop("Invalid perturbation targets.")
annotations <- fread(annotation_path)
annotations[, method := fifelse(method == "SCENIC+", "SCENICplus", method)]
annotations <- annotations[method %chin% methods]

network_cache <- list()
for (method_name in methods) {
  network_cache[[method_name]] <- list()
  for (day_slug in day_slugs) {
    path <- file.path(network_root, method_name, day_slug, "tf_gene.tsv.gz")
    network_cache[[method_name]][[day_slug]] <- collapse_tf_gene(path)
  }
}

all_members <- list()
all_enrichment <- list()
all_neighborhood_metrics <- list()
all_atac_members <- list()
counter <- 0L

for (day_slug in day_slugs) {
  day_label <- day_labels[[day_slug]]
  for (target_tf in targets) {
    rna_path <- file.path(truth_root, "rna_effects", paste0(day_slug, "__", target_tf, ".tsv.gz"))
    atac_path <- file.path(truth_root, "atac_effects", paste0(day_slug, "__", target_tf, ".tsv.gz"))
    rna_truth <- fread(rna_path)
    atac_truth <- fread(atac_path)
    rna_stats <- perturbation_response_stat(rna_truth)
    atac_stats <- perturbation_response_stat(atac_truth)

    pathways <- list()
    pathway_methods <- character()
    pathway_orders <- integer()
    for (method_name in methods) {
      neighborhood <- build_neighborhood(network_cache[[method_name]][[day_slug]], target_tf)
      if (nrow(neighborhood)) {
        neighborhood[, `:=`(method = method_name, perturbation_regulator = target_tf, timepoint = day_label)]
        all_members[[length(all_members) + 1L]] <- neighborhood
        for (order_value in 1:2) {
          member_set <- neighborhood[order == order_value]
          if (!nrow(member_set)) next
          key <- paste(method_name, "RNA", paste0("order", order_value), sep = "|")
          pathways[[key]] <- unique(member_set$target)
          pathway_methods[[key]] <- method_name
          pathway_orders[[key]] <- order_value
          metric <- evaluate_set(rna_truth, member_set)
          if (nrow(metric)) {
            metric[, `:=`(
              method = method_name,
              perturbation_regulator = target_tf,
              timepoint = day_label,
              modality = "RNA",
              neighborhood_order = order_value
            )]
            all_neighborhood_metrics[[length(all_neighborhood_metrics) + 1L]] <- metric
          }
        }
      }
    }
    if (length(pathways)) {
      enr <- run_enrichment(rna_stats, pathways, "RNA", pathway_methods, pathway_orders)
      if (nrow(enr)) {
        enr[, `:=`(perturbation_regulator = target_tf, timepoint = day_label)]
        all_enrichment[[length(all_enrichment) + 1L]] <- enr
      }
    }

    atac_pathways <- list()
    atac_methods <- character()
    atac_orders <- integer()
    for (method_name in methods) {
      mapped <- unique(annotations[
        method == method_name & perturbation_regulator == target_tf & timepoint == day_label &
          atac_evaluable == TRUE & !is.na(atac_truth_feature),
        .(target = atac_truth_feature)
      ])
      if (!nrow(mapped)) next
      mapped[, `:=`(
        order = 1L,
        intermediate = NA_character_,
        path_direction = NA_integer_,
        paths = 1L,
        method = method_name,
        perturbation_regulator = target_tf,
        timepoint = day_label
      )]
      all_atac_members[[length(all_atac_members) + 1L]] <- mapped
      key <- paste(method_name, "ATAC", "order1", sep = "|")
      atac_pathways[[key]] <- mapped$target
      atac_methods[[key]] <- method_name
      atac_orders[[key]] <- 1L
      metric <- evaluate_set(atac_truth, mapped)
      if (nrow(metric)) {
        metric[, `:=`(
          method = method_name,
          perturbation_regulator = target_tf,
          timepoint = day_label,
          modality = "ATAC",
          neighborhood_order = 1L
        )]
        all_neighborhood_metrics[[length(all_neighborhood_metrics) + 1L]] <- metric
      }
    }
    if (length(atac_pathways)) {
      enr <- run_enrichment(atac_stats, atac_pathways, "ATAC", atac_methods, atac_orders)
      if (nrow(enr)) {
        enr[, `:=`(perturbation_regulator = target_tf, timepoint = day_label)]
        all_enrichment[[length(all_enrichment) + 1L]] <- enr
      }
    }
    counter <- counter + 1L
    if (counter %% 10L == 0L) message("Completed ", counter, "/", length(targets) * length(day_slugs), " perturbation strata")
  }
}

members <- rbindlist(all_members, fill = TRUE)
atac_members <- rbindlist(all_atac_members, fill = TRUE)
enrichment <- rbindlist(all_enrichment, fill = TRUE)
metrics <- rbindlist(all_neighborhood_metrics, fill = TRUE)
if (!nrow(enrichment) || !nrow(metrics)) stop("No perturbation neighborhood results were produced.")
enrichment[, padj_global := p.adjust(pval, method = "BH"), by = .(modality, neighborhood_order)]
metrics[, response_p_adjust := p.adjust(response_p_value, method = "BH"), by = .(modality, neighborhood_order)]

direct_enrichment <- enrichment[neighborhood_order == 1L, .(
  method, perturbation_regulator, timepoint, modality,
  set_size = size, enrichment_score = ES, normalized_enrichment_score = NES,
  enrichment_p_value = pval, enrichment_p_adjust = padj_global, leading_edge
)]
joint <- dcast(
  direct_enrichment,
  method + perturbation_regulator + timepoint ~ modality,
  value.var = "normalized_enrichment_score"
)
if (!all(c("RNA", "ATAC") %in% names(joint))) stop("Joint enrichment modalities are incomplete.")
joint[, joint_NES := rowMeans(.SD, na.rm = FALSE), .SDcols = c("RNA", "ATAC")]

tf_summary <- joint[, .(
  timepoints = .N,
  median_RNA_NES = median(RNA),
  median_ATAC_NES = median(ATAC),
  median_joint_NES = median(joint_NES),
  recurrent_positive_timepoints = sum(joint_NES > 0)
), by = .(method, perturbation_regulator)]

n_methods <- length(methods)
common_tf <- tf_summary[, .(methods = uniqueN(method), min_timepoints = min(timepoints)), by = perturbation_regulator][methods == n_methods & min_timepoints >= 3L, perturbation_regulator]
hero_candidates <- metrics[
  method == "multiCSN" & modality == "RNA" & neighborhood_order %in% 1:2 &
    perturbation_regulator %chin% common_tf,
  .(
    responsive = sum(responsive),
    direction_consistent = sum(direction_consistent),
    evaluated_timepoints = uniqueN(timepoint),
    median_odds_ratio = median(response_odds_ratio, na.rm = TRUE)
  ),
  by = perturbation_regulator
]
setorder(hero_candidates, -direction_consistent, -responsive, -median_odds_ratio, perturbation_regulator)
if (!nrow(hero_candidates)) stop("No common hero TF candidate.")
hero_tf <- hero_candidates$perturbation_regulator[[1L]]
hero_day <- joint[method == "multiCSN" & perturbation_regulator == hero_tf][order(-joint_NES, timepoint), timepoint][1L]

hero_members <- members[method == "multiCSN" & perturbation_regulator == hero_tf & timepoint == hero_day]
hero_truth_path <- file.path(truth_root, "rna_effects", paste0(names(day_labels)[match(hero_day, day_labels)], "__", hero_tf, ".tsv.gz"))
hero_truth <- fread(hero_truth_path)
hero_paths <- merge(hero_members, hero_truth[, .(target = feature, rna_logFC = logFC, rna_FDR = FDR, rna_PValue = PValue)], by = "target")
hero_paths <- hero_paths[rna_FDR < 0.05 & is.finite(path_direction) & sign(rna_logFC) == -path_direction]
hero_paths[, response_rank := frank(rna_PValue, ties.method = "first")]
setorder(hero_paths, order, response_rank, target)
hero_paths <- hero_paths[, head(.SD, if (order[1L] == 1L) 6L else 12L), by = order]
hero_paths[, selected_intermediate := fifelse(order == 2L, tstrsplit(intermediate, ";", fixed = TRUE, keep = 1L)[[1L]], NA_character_)]

hero_curve <- list()
hero_day_slug <- names(day_labels)[match(hero_day, day_labels)]
hero_rna <- fread(file.path(truth_root, "rna_effects", paste0(hero_day_slug, "__", hero_tf, ".tsv.gz")))
hero_atac <- fread(file.path(truth_root, "atac_effects", paste0(hero_day_slug, "__", hero_tf, ".tsv.gz")))
for (modality in c("RNA", "ATAC")) {
  truth <- if (modality == "RNA") hero_rna else hero_atac
  stats <- perturbation_response_stat(truth)
  for (method_name in methods) {
    set_features <- if (modality == "RNA") {
      members[method == method_name & perturbation_regulator == hero_tf & timepoint == hero_day & order == 1L, unique(target)]
    } else {
      atac_members[method == method_name & perturbation_regulator == hero_tf & timepoint == hero_day, unique(target)]
    }
    hits <- names(stats) %chin% set_features
    if (sum(hits) < 5L) next
    hit_weight <- abs(stats) * hits
    miss_penalty <- (!hits) / sum(!hits)
    hit_reward <- if (sum(hit_weight) > 0) hit_weight / sum(hit_weight) else hits / sum(hits)
    running <- cumsum(hit_reward - miss_penalty)
    hero_curve[[length(hero_curve) + 1L]] <- data.table(
      method = method_name,
      modality = modality,
      perturbation_regulator = hero_tf,
      timepoint = hero_day,
      rank = seq_along(running),
      rank_fraction = seq_along(running) / length(running),
      running_enrichment = running,
      hit = hits,
      feature = names(stats)
    )
  }
}
hero_curve <- rbindlist(hero_curve)

fwrite(members, file.path(stage, "rna_network_neighborhood_members.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(atac_members, file.path(stage, "atac_direct_region_members.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(enrichment, file.path(stage, "full_background_enrichment.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(metrics, file.path(stage, "neighborhood_response_metrics.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")
fwrite(joint, file.path(stage, "joint_direct_enrichment.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(tf_summary, file.path(stage, "TF_enrichment_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(hero_candidates, file.path(stage, "hero_TF_selection.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(hero_paths, file.path(stage, "hero_direction_consistent_paths.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(hero_curve, file.path(stage, "hero_full_background_enrichment_curves.tsv.gz"), sep = "\t", quote = FALSE, na = "NA")

publish_result_files(stage, output)
