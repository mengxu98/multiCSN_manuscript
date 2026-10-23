#!/usr/bin/env Rscript

source("code/functions/ranking_metrics.R")

source("code/functions/table_io.R")

suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)

evaluate_linger_multilayer <- function(mode = "binding", networks = "runs/gse274113/LINGER",
  endpoints = "results/gse274113/networks", truth = "data/gse274113/heldout_pseudobulk_full",
  evidence = "data/gse274113/external_resources", object = "data/gse274113/objects/GSE274113_control_full",
  chain_helper = "code/gse274113/evaluate_chain_response.R", output = NULL) {
  mode <- match.arg(mode, c("binding", "contact", "chain"))
  if (is.null(output)) output <- "results/gse274113/evidence"
  arguments <- list(mode = mode, networks = networks, endpoints = endpoints, truth = truth,
                    evidence = evidence, object = object, chain_helper = chain_helper, output = output)
network_root <- normalizePath(arguments$networks, mustWork = TRUE)
truth_root <- normalizePath(arguments$truth, mustWork = TRUE)
endpoint_root <- normalizePath(arguments$endpoints, mustWork = FALSE)
output_root <- normalizePath(arguments$output, mustWork = FALSE)
timepoints <- c("day 7", "day 9", "day 11", "day 14")
effect_summary_path <- file.path(truth_root, "effect_model_summary.tsv")
effect_summary <- fread(effect_summary_path)
targets <- sort(unique(as.character(effect_summary$target)))
if (length(targets) != 19L ||
    !setequal(unique(effect_summary$timepoint), timepoints)) {
  stop("Held-out truth failed the 19-TF four-timepoint gate.", call. = FALSE)
}

prepare_endpoints <- function() {
  required <- unlist(lapply(gsub(" ", "_", timepoints, fixed = TRUE), function(slug) {
    file.path(endpoint_root, "LINGER", slug, c("tf_region.tsv.gz", "region_gene.tsv.gz"))
  }))
  linger_root <- file.path(endpoint_root, "LINGER")
  if (all(file.exists(required))) return(invisible(NULL))
  dir.create(endpoint_root, recursive = TRUE, showWarnings = FALSE)
  stage <- paste0(linger_root, ".staging_", Sys.getpid())
  if (file.exists(stage)) stop("Staging output already exists: ", stage)
  dir.create(stage, recursive = FALSE)
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)
  summary_rows <- list()
  for (timepoint in timepoints) {
    slug <- gsub(" ", "_", timepoint, fixed = TRUE)
    source_dir <- file.path(network_root, slug, "networks")
    destination <- file.path(stage, slug)
    dir.create(destination, recursive = TRUE)
    tf_region_path <- file.path(source_dir, "tf_region.tsv.gz")
    region_gene_path <- file.path(source_dir, "region_gene.tsv.gz")
    summary_path <- file.path(source_dir, "network_summary.tsv")
    if (!all(file.exists(c(tf_region_path, region_gene_path, summary_path)))) {
      stop("Missing LINGER endpoint at ", slug, call. = FALSE)
    }
    tf_region <- fread(tf_region_path)
    if (!all(c("region", "regulator", "score") %in% names(tf_region)) ||
        any(!is.finite(tf_region$score)) || any(tf_region$score == 0)) {
      stop("Invalid LINGER TF-region endpoint at ", slug, call. = FALSE)
    }
    raw_rows <- nrow(tf_region)
    tf_region <- tf_region[regulator %chin% targets]
    evaluated_rows_before_deduplication <- nrow(tf_region)
    tf_region <- tf_region[, .(score = max(score)), by = .(region, regulator)]
    duplicate_tf_region_rows <- evaluated_rows_before_deduplication - nrow(tf_region)
    write_tsv(tf_region, file.path(destination, "tf_region.tsv.gz"))
    region_gene <- fread(region_gene_path)
    if (!all(c("region", "target", "score") %in% names(region_gene)) ||
        any(!is.finite(region_gene$score)) || any(region_gene$score == 0)) {
      stop("Invalid LINGER region-gene endpoint at ", slug, call. = FALSE)
    }
    region_gene_rows_before_deduplication <- nrow(region_gene)
    region_gene <- region_gene[, .(score = max(score)), by = .(region, target)]
    duplicate_region_gene_rows <- region_gene_rows_before_deduplication - nrow(region_gene)
    write_tsv(region_gene, file.path(destination, "region_gene.tsv.gz"))
    endpoint_summary <- fread(summary_path)
    summary_rows[[slug]] <- data.table(
      timepoint = timepoint,
      raw_tf_region_relations = raw_rows,
      evaluated_tf_region_relations = nrow(tf_region),
      duplicate_tf_region_rows_collapsed = duplicate_tf_region_rows,
      perturbed_tfs_available = uniqueN(tf_region$regulator),
      missing_perturbed_tfs = paste(setdiff(targets, unique(tf_region$regulator)), collapse = ","),
      region_gene_relations = nrow(region_gene),
      duplicate_region_gene_rows_collapsed = duplicate_region_gene_rows,
      selection_policy = "all_finite_nonzero_LINGER_native_scores"
    )
    if (endpoint_summary[endpoint == "region_gene", relations] !=
        region_gene_rows_before_deduplication) {
      stop("LINGER region-gene row accounting changed at ", slug, call. = FALSE)
    }
    rm(tf_region, region_gene)
    gc(verbose = FALSE)
  }
  write_tsv(rbindlist(summary_rows), file.path(stage, "endpoint_summary.tsv"))
  dir.create(linger_root, recursive = TRUE, showWarnings = FALSE)
  for (relative in list.files(stage, recursive = TRUE)) {
    destination <- file.path(linger_root, relative)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!file.rename(file.path(stage, relative), destination)) {
      stop("Unable to publish prepared endpoint: ", destination)
    }
  }
}

run_endpoint_evaluator_for_linger <- function(mode) {
  evaluator_name <- switch(mode,
    binding = "evaluate_tf_region_binding",
    contact = "evaluate_region_gene_contact",
    stop("Unknown endpoint evaluation mode: ", mode)
  )
  environment <- new.env(parent = globalenv())
  sys.source(file.path("code/gse274113", paste0(evaluator_name, ".R")), envir = environment)
  settings <- list(
    endpoints = endpoint_root,
    evidence = normalizePath(arguments$evidence, mustWork = TRUE),
    output = output_root,
    methods = "LINGER"
  )
  if (mode == "binding") {
    settings$truth <- truth_root
  } else {
    settings$object <- normalizePath(arguments$object, mustWork = TRUE)
  }
  environment[[evaluator_name]](settings)
}

run_chain_evaluation <- function() {
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(IRanges)
  })
  source(normalizePath(arguments$chain_helper, mustWork = TRUE), local = TRUE)
  dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(".result-files-", tmpdir = dirname(output_root))
  dir.create(stage, recursive = FALSE)
  metric_rows <- list()
  cursor <- 1L
  for (timepoint in timepoints) {
    slug <- gsub(" ", "_", timepoint, fixed = TRUE)
    tf_region_path <- file.path(endpoint_root, "LINGER", slug, "tf_region.tsv.gz")
    region_gene_path <- file.path(endpoint_root, "LINGER", slug, "region_gene.tsv.gz")
    tf_region <- fread(tf_region_path)[, .(
      regulator = as.character(regulator), region = as.character(region)
    )]
    region_gene <- fread(region_gene_path)[, .(
      region = as.character(region), target = as.character(target)
    )]
    if (anyDuplicated(tf_region) || anyDuplicated(region_gene)) {
      stop("Duplicated LINGER endpoint keys at ", slug, call. = FALSE)
    }
    setkey(region_gene, region)
    evaluated_targets <- parallel::mclapply(targets, function(regulator) {
      regulator_id <- regulator
      tf_regions <- tf_region[regulator == regulator_id, .(region)]
      if (nrow(tf_regions)) {
        chains <- region_gene[tf_regions, on = "region", nomatch = 0L, allow.cartesian = TRUE][
          , .(
            regulator = regulator_id, region, target,
            chain_direction = NA_real_, tf_region_direction = NA_real_,
            region_gene_direction = NA_real_
          )
        ]
      } else {
        chains <- data.table(
          regulator = character(), region = character(), target = character(),
          chain_direction = numeric(), tf_region_direction = numeric(),
          region_gene_direction = numeric()
        )
      }
      truth_slug <- paste(slug, regulator, sep = "__")
      atac_path <- file.path(truth_root, "atac_effects", paste0(truth_slug, ".tsv.gz"))
      rna_path <- file.path(truth_root, "rna_effects", paste0(truth_slug, ".tsv.gz"))
      evaluated <- evaluate_selected_chain_response(
        chains, fread(atac_path), fread(rna_path),
        "LINGER", regulator, timepoint
      )$metrics
      list(metrics = evaluated, triplet_count = nrow(chains))
    }, mc.cores = min(8L, length(targets)), mc.preschedule = FALSE)
    if (any(vapply(evaluated_targets, inherits, logical(1), "try-error"))) {
      stop("Parallel LINGER chain evaluation failed at ", timepoint)
    }
    triplet_count <- sum(vapply(
      evaluated_targets, `[[`, numeric(1), "triplet_count"
    ))
    for (evaluated in evaluated_targets) {
      metric_rows[[cursor]] <- evaluated$metrics
      cursor <- cursor + 1L
    }
    cat("chain_timepoint\tPASS\t", timepoint, "\t", triplet_count, "\n", sep = "")
    rm(tf_region, region_gene, evaluated_targets)
    gc(verbose = FALSE)
  }
  strata <- rbindlist(metric_rows, use.names = TRUE, fill = TRUE)
  if (nrow(strata) != 76L || anyDuplicated(strata[, .(method, regulator, timepoint)])) {
    stop("LINGER chain response failed the 19-TF four-timepoint gate.")
  }
  clusters <- summarize_chain_responses(strata)
  write_tsv(strata, file.path(stage, "chain_response_stratum_metrics.tsv.gz"))
  write_tsv(clusters, file.path(stage, "chain_response_TF_cluster_metrics.tsv"))
  publish_result_files(stage, output_root, prefix = "linger_")
  cat("linger_chain_evaluation\tPASS\t", output_root, "\n", sep = "")
}

prepare_endpoints()
if (arguments$mode %chin% c("binding", "contact")) {
  run_endpoint_evaluator_for_linger(arguments$mode)
} else {
  run_chain_evaluation()
}

}

run_linger_binding <- function() evaluate_linger_multilayer("binding")
run_linger_contact <- function() evaluate_linger_multilayer("contact")
run_linger_chain <- function() evaluate_linger_multilayer("chain")

if (sys.nframe() == 0L) run_linger_binding()
