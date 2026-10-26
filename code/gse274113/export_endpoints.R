#!/usr/bin/env Rscript
# Export multiCSN, Pando and SCENIC+ endpoints with a shared score schema.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
source("code/functions/network_io.R")
require_columns <- function(table, columns, label) {
  if (!all(columns %in% names(table))) {
    stop(label, " schema changed; missing: ",
         paste(setdiff(columns, names(table)), collapse = ", "))
  }
}

validate_standard <- function(table, keys, label) {
  required <- c(keys, "rank_magnitude", "direction_sign", "score_available")
  require_columns(table, required, label)
  if (anyDuplicated(table[, ..keys]) ||
      any(table$score_available &
          (!is.finite(table$rank_magnitude) | table$rank_magnitude < 0)) ||
      any(!table$score_available & table$rank_magnitude != 0) ||
      any(is.finite(table$direction_sign) &
          !table$direction_sign %in% c(-1, 0, 1))) {
    stop("Invalid standardized endpoint: ", label)
  }
  invisible(TRUE)
}

collapse_best <- function(table, keys, score_column, direction_column, label) {
  table <- as.data.table(table)
  table[, method_score := as.numeric(get(score_column))]
  table[, method_direction := sign(as.numeric(get(direction_column)))]
  output <- table[, {
    available <- is.finite(method_score) & method_score >= 0
    if (!any(available)) {
      list(rank_magnitude = 0, direction_sign = NA_real_, score_available = FALSE)
    } else {
      best <- max(method_score[available])
      top_sign <- unique(method_direction[available & method_score == best])
      top_sign <- top_sign[is.finite(top_sign)]
      if (length(top_sign) > 1L) {
        stop(label, " has conflicting directions among top-score duplicates.")
      }
      list(
        rank_magnitude = best,
        direction_sign = if (length(top_sign)) top_sign[[1]] else NA_real_,
        score_available = TRUE
      )
    }
  }, by = keys]
  validate_standard(output, keys, label)
  output
}

pvalue_score <- function(p) {
  p <- as.numeric(p)
  output <- numeric(length(p))
  keep <- is.finite(p) & p >= 0 & p <= 1
  output[keep] <- -log10(pmax(p[keep], .Machine$double.xmin))
  list(score = output, available = keep)
}


export_endpoints <- function(method, input, output) {
  if (!method %in% c("multiCSN", "Pando", "SCENICplus")) stop("Unknown method.")
  input <- normalizePath(input, mustWork = TRUE)
  output <- normalizePath(output, mustWork = FALSE)
  if (file.exists(output)) stop("Output already exists: ", output)
  stage_root <- paste0(output, ".staging")
  if (file.exists(stage_root)) stop("Inspect existing staging output before retrying: ", stage_root)
  dir.create(stage_root, recursive = TRUE)
  write_endpoint <- function(method, slug, endpoint, table, keys, score_semantics) {
    validate_standard(table, keys, paste(method, slug, endpoint))
    write_tsv(table, file.path(stage_root, slug, paste0(gsub("-", "_", tolower(endpoint)), ".tsv.gz")))
  }
  for (slug in c("day_7", "day_9", "day_11", "day_14")) {
    dir.create(file.path(stage_root, slug))
    roots <- list(multicsn = file.path(input, slug), pando = input, scenicplus = input)
    if (method == "multiCSN") {
  multi_dir <- roots$multicsn
  multi_paths <- file.path(multi_dir, c("tf_gene.tsv", "tf_region.tsv", "region_gene.tsv", "triplets.tsv"))
  if (!all(file.exists(multi_paths))) stop("Missing multiCSN endpoints: ", multi_dir)
  specs <- list(tf_gene = c("regulator", "target"), tf_region = c("regulator", "region"),
                region_gene = c("region", "target"), triplet = c("regulator", "region", "target"))
  for (i in seq_along(specs)) {
    table <- read_fit_checkpoint(multi_paths[[i]])
    keys <- specs[[i]]
    require_columns(table, c(keys, "weight"), names(specs)[[i]])
    standardized <- table[, c(setNames(lapply(.SD, as.character), keys),
      list(rank_magnitude = abs(as.numeric(weight)), direction_sign = sign(as.numeric(weight)),
           score_available = TRUE)), .SDcols = keys]
    write_endpoint("multiCSN", slug, names(specs)[[i]], standardized, keys, "within-timepoint ordinal weight")
  }
    } else if (method == "Pando") {
  pando_dir <- file.path(roots$pando, slug)
  pando_paths <- file.path(pando_dir, c(
    "native_module_metadata.tsv", "retained_native_triplets.tsv.gz",
    "retained_tf_region_binary.tsv.gz", "retained_region_gene_binary.tsv.gz"
  ))
  if (!all(file.exists(pando_paths))) stop("Missing Pando endpoints at ", slug)
  pando_tg_raw <- read_fit_checkpoint(pando_paths[[1]])
  require_columns(pando_tg_raw, c("tf", "target", "estimate", "padj"),
                  "Pando TF-gene")
  pando_score <- pvalue_score(pando_tg_raw$padj)
  pando_tg <- pando_tg_raw[, .(regulator = tf, target)]
  pando_tg[, `:=`(
    rank_magnitude = pando_score$score,
    direction_sign = sign(as.numeric(pando_tg_raw$estimate)),
    score_available = pando_score$available
  )]
  pando_triplet_raw <- read_fit_checkpoint(pando_paths[[2]])
  require_columns(pando_triplet_raw,
                  c("tf", "region", "target", "estimate", "padj"),
                  "Pando triplet")
  pando_triplet_score <- pvalue_score(pando_triplet_raw$padj)
  pando_triplet <- pando_triplet_raw[, .(regulator = tf, region, target)]
  pando_triplet[, `:=`(
    rank_magnitude = pando_triplet_score$score,
    direction_sign = sign(as.numeric(pando_triplet_raw$estimate)),
    score_available = pando_triplet_score$available
  )]
  pando_tr <- read_fit_checkpoint(pando_paths[[3]])
  require_columns(pando_tr, c("tf", "region"), "Pando TF-region")
  pando_tr <- pando_tr[, .(
    regulator = tf, region, rank_magnitude = 0,
    direction_sign = NA_real_, score_available = FALSE
  )]
  pando_rg <- read_fit_checkpoint(pando_paths[[4]])
  require_columns(pando_rg, c("region", "target"), "Pando region-gene")
  pando_rg <- pando_rg[, .(
    region, target, rank_magnitude = 0,
    direction_sign = NA_real_, score_available = FALSE
  )]
  write_endpoint("Pando", slug, "TF-gene", pando_tg,
                 c("regulator", "target"), "adjusted-p-value ranking")
  write_endpoint("Pando", slug, "TF-region", pando_tr,
                 c("regulator", "region"), "binary projection of retained triplet terms")
  write_endpoint("Pando", slug, "region-gene", pando_rg,
                 c("region", "target"), "binary projection of retained triplet terms")
  write_endpoint("Pando", slug, "triplet", pando_triplet,
                 c("regulator", "region", "target"),
                 "interaction-term adjusted-p-value ranking")
    } else {
  scenic_dir <- file.path(roots$scenicplus, slug)
  scenic_path <- file.path(scenic_dir, "eRegulon_direct.tsv")
  if (!file.exists(scenic_path)) {
    stop("Missing SCENIC+ direct eRegulons at ", slug)
  }
  scenic <- read_fit_checkpoint(scenic_path)
  scenic_required <- c(
    "TF", "Region", "Gene", "importance_R2G", "rho_R2G",
    "importance_TF2G", "rho_TF2G", "triplet_rank"
  )
  require_columns(scenic, scenic_required, "SCENIC+ direct eRegulon")
  scenic[, `:=`(
    importance_R2G = as.numeric(importance_R2G),
    rho_R2G = as.numeric(rho_R2G),
    importance_TF2G = as.numeric(importance_TF2G),
    rho_TF2G = as.numeric(rho_TF2G),
    triplet_rank = as.numeric(triplet_rank)
  )]
  if (any(!is.finite(scenic$triplet_rank) | scenic$triplet_rank < 0)) {
    stop("SCENIC+ triplet_rank is invalid at ", slug)
  }
  scenic_tg <- collapse_best(
    scenic, c("TF", "Gene"), "importance_TF2G", "rho_TF2G",
    paste("SCENIC+ TF-gene", slug)
  )
  setnames(scenic_tg, c("TF", "Gene"), c("regulator", "target"))
  scenic_rg <- collapse_best(
    scenic, c("Region", "Gene"), "importance_R2G", "rho_R2G",
    paste("SCENIC+ region-gene", slug)
  )
  setnames(scenic_rg, c("Region", "Gene"), c("region", "target"))
  scenic_tr <- unique(scenic[, .(regulator = TF, region = Region)])
  scenic_tr[, `:=`(
    rank_magnitude = 0, direction_sign = NA_real_, score_available = FALSE
  )]
  scenic[, triplet_magnitude := max(triplet_rank) + 1 - triplet_rank]
  scenic_triplet <- collapse_best(
    scenic, c("TF", "Region", "Gene"), "triplet_magnitude", "rho_TF2G",
    paste("SCENIC+ triplet", slug)
  )
  setnames(
    scenic_triplet,
    c("TF", "Region", "Gene"),
    c("regulator", "region", "target")
  )
  write_endpoint("SCENICplus", slug, "TF-gene", scenic_tg,
                 c("regulator", "target"), "TF-to-gene importance")
  write_endpoint("SCENICplus", slug, "TF-region", scenic_tr,
                 c("regulator", "region"), "binary direct-cistrome membership")
  write_endpoint("SCENICplus", slug, "region-gene", scenic_rg,
                 c("region", "target"), "region-to-gene importance")
  write_endpoint("SCENICplus", slug, "triplet", scenic_triplet,
                 c("regulator", "region", "target"), "triplet_rank")
    }
  }
  if (!file.rename(stage_root, output)) stop("Unable to promote canonical endpoints.")
}
if (sys.nframe() == 0L) {
  # Default: export the three formal methods from their method-specific fit trees.
  for (method in c("multiCSN", "Pando", "SCENICplus")) {
    export_endpoints(
      method,
      file.path("runs", "gse274113", "method_fits", method),
      file.path("results", "gse274113", "networks", method)
    )
  }
}
