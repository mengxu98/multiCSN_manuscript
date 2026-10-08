#!/usr/bin/env Rscript

source("code/functions/table_io.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
})
data.table::setDTthreads(1L)


regions_to_table <- function(values) {
  values <- as.character(values)
  matched <- regexec("^(.+?)[-:]([0-9]+)-([0-9]+)$", values, perl = TRUE)
  pieces <- regmatches(values, matched)
  valid <- lengths(pieces) == 4L
  if (!all(valid)) stop("Unparseable region coordinates.")
  data.table(
    region = values,
    chrom = vapply(pieces, `[[`, character(1), 2L),
    start = as.integer(vapply(pieces, `[[`, character(1), 3L)),
    end = as.integer(vapply(pieces, `[[`, character(1), 4L))
  )
}

table_to_granges <- function(table, chrom = "chrom", start = "start", end = "end") {
  GenomicRanges::GRanges(
    as.character(table[[chrom]]),
    IRanges::IRanges(start = as.integer(table[[start]]), end = as.integer(table[[end]]))
  )
}

distance_breaks <- c(
  -Inf, 10000, 25000, 50000, 100000, 250000,
  500000, 1000000, Inf
)
distance_labels <- c(
  "0-10kb", "10-25kb", "25-50kb", "50-100kb",
  "100-250kb", "250-500kb", "500kb-1Mb", ">1Mb"
)
distance_bin <- function(distance) {
  as.character(cut(
    distance,
    breaks = distance_breaks, labels = distance_labels,
    right = FALSE, include.lowest = TRUE
  ))
}

# ---- Parameters ----
# Prepared HSPC evidence was rebuilt on 2026-09-14 with recorded UCSC chains.
# See data/gse274113/raw_external_downloads/README.md for provenance.
# object: NT control directory; output: region_gene_HSPC_contact.
evaluate_region_gene_contact <- function(arguments = list(
  endpoints = "results/gse274113/networks",
  evidence = "data/gse274113/external_resources",
  object = "data/gse274113/objects/GSE274113_control_full",
  output = "results/gse274113/multiomic_endpoint_evidence/region_gene_HSPC_contact",
  methods = NULL
)) {
endpoint_root <- normalizePath(arguments$endpoints, mustWork = TRUE)
evidence_root <- normalizePath(arguments$evidence, mustWork = TRUE)
object_root <- normalizePath(arguments$object, mustWork = TRUE)
output_root <- normalizePath(arguments$output, mustWork = FALSE)
if (file.exists(output_root)) {
  stop("Output already exists and will not be overwritten: ", output_root)
}
dir.create(dirname(output_root), recursive = TRUE, showWarnings = FALSE)
stage_root <- paste0(
  output_root, ".staging_", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
)
dir.create(stage_root, recursive = FALSE)

loop_path <- file.path(evidence_root, "HSPC_loops_hg38.tsv.gz")
if (!file.exists(loop_path)) stop("Prepared HSPC loops are required.")
loops <- fread(loop_path)
loop_columns <- c(
  "loop_id", "anchor1_chrom", "anchor1_start", "anchor1_end",
  "anchor2_chrom", "anchor2_start", "anchor2_end"
)
if (!all(loop_columns %in% names(loops))) stop("HSPC loop schema changed.")
anchor1 <- table_to_granges(
  loops, "anchor1_chrom", "anchor1_start", "anchor1_end"
)
anchor2 <- table_to_granges(
  loops, "anchor2_chrom", "anchor2_start", "anchor2_end"
)

genes <- ensembldb::genes(
  EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
  return.type = "GRanges"
)
seqlevel_map <- stats::setNames(
  c(paste0("chr", 1:22), "chrX", "chrY", "chrM"),
  c(as.character(1:22), "X", "Y", "MT")
)
seqlevel_map <- seqlevel_map[names(seqlevel_map) %in% GenomeInfoDb::seqlevels(genes)]
genes <- GenomeInfoDb::renameSeqlevels(genes, seqlevel_map)
gene_name <- as.character(S4Vectors::mcols(genes)$gene_name)
keep_gene <- !is.na(gene_name) & nzchar(gene_name) &
  as.character(GenomeInfoDb::seqnames(genes)) %in% c(
    paste0("chr", 1:22), "chrX", "chrY"
  )
genes <- genes[keep_gene]
gene_name <- gene_name[keep_gene]
tss_position <- ifelse(as.character(GenomicRanges::strand(genes)) == "-",
  GenomicRanges::end(genes), GenomicRanges::start(genes)
)
tss <- GenomicRanges::GRanges(
  GenomeInfoDb::seqnames(genes), IRanges::IRanges(tss_position, width = 1)
)
S4Vectors::mcols(tss)$gene_name <- gene_name
tss_table <- data.table(
  target = gene_name,
  chrom = as.character(GenomeInfoDb::seqnames(tss)),
  tss = tss_position
)
tss_table <- unique(tss_table)

support_rows <- list()
support_index <- 0L
hits1 <- GenomicRanges::findOverlaps(tss, anchor1, ignore.strand = TRUE)
if (length(hits1)) {
  support_index <- support_index + 1L
  query <- S4Vectors::queryHits(hits1)
  subject <- S4Vectors::subjectHits(hits1)
  support_rows[[support_index]] <- data.table(
    target = gene_name[query],
    chrom = as.character(GenomeInfoDb::seqnames(anchor2[subject])),
    start = GenomicRanges::start(anchor2[subject]),
    end = GenomicRanges::end(anchor2[subject]),
    loop_id = loops$loop_id[subject]
  )
}
hits2 <- GenomicRanges::findOverlaps(tss, anchor2, ignore.strand = TRUE)
if (length(hits2)) {
  support_index <- support_index + 1L
  query <- S4Vectors::queryHits(hits2)
  subject <- S4Vectors::subjectHits(hits2)
  support_rows[[support_index]] <- data.table(
    target = gene_name[query],
    chrom = as.character(GenomeInfoDb::seqnames(anchor1[subject])),
    start = GenomicRanges::start(anchor1[subject]),
    end = GenomicRanges::end(anchor1[subject]),
    loop_id = loops$loop_id[subject]
  )
}
support_table <- unique(rbindlist(support_rows, use.names = TRUE))
if (!nrow(support_table)) stop("No annotated TSS maps to an HSPC loop anchor.")
support_ranges <- lapply(split(support_table, support_table$target), function(table) {
  GenomicRanges::reduce(table_to_granges(table))
})

methods <- if (is.null(arguments$methods)) {
  c("multiCSN", "SCENICplus", "Pando")
} else {
  strsplit(arguments$methods, ",", fixed = TRUE)[[1]]
}
if (!length(methods) || anyDuplicated(methods) || any(!nzchar(methods))) {
  stop("methods must contain unique method names.")
}
timepoints <- c("day 7", "day 9", "day 11", "day 14")
result_rows <- list()
edge_rows <- list()
result_index <- 0L
edge_index <- 0L

for (timepoint in timepoints) {
  slug <- gsub(" ", "_", timepoint, fixed = TRUE)
  accessible_path <- file.path(
    object_root, paste0("accessible_consensus_peaks_", slug, ".tsv")
  )
  if (!file.exists(accessible_path)) stop("Missing accessible peak universe: ", accessible_path)
  accessible_table <- regions_to_table(fread(accessible_path)$peak)
  accessible_ranges <- table_to_granges(accessible_table)
  accessible_midpoint <- floor((accessible_table$start + accessible_table$end) / 2)

  prediction_tables <- list()
  target_union <- character()
  for (method in methods) {
    path <- file.path(endpoint_root, method, slug, "region_gene.tsv.gz")
    if (!file.exists(path)) stop("Missing region-gene endpoint: ", path)
    prediction <- fread(path)
    if (!all(c("region", "target") %in% names(prediction)) ||
      anyDuplicated(prediction[, .(region, target)])) {
      stop("Invalid region-gene endpoint: ", method, " ", timepoint)
    }
    prediction_tables[[method]] <- prediction
    target_union <- union(target_union, prediction$target)
  }
  target_union <- intersect(target_union, names(support_ranges))
  background_cache <- new.env(parent = emptyenv())

  build_background <- function(target) {
    target_id <- target
    if (exists(target, envir = background_cache, inherits = FALSE)) {
      return(get(target, envir = background_cache, inherits = FALSE))
    }
    target_tss <- tss_table[target == target_id]
    keep <- accessible_table$chrom %chin% target_tss$chrom
    peaks <- accessible_ranges[keep]
    table <- accessible_table[keep]
    midpoint <- accessible_midpoint[keep]
    distance <- rep(Inf, nrow(table))
    for (chromosome in unique(table$chrom)) {
      peak_index <- which(table$chrom == chromosome)
      positions <- target_tss[chrom == chromosome, tss]
      if (!length(positions)) next
      distance[peak_index] <- vapply(midpoint[peak_index], function(value) {
        min(abs(value - positions))
      }, numeric(1))
    }
    supported <- GenomicRanges::countOverlaps(
      peaks, support_ranges[[target]],
      ignore.strand = TRUE
    ) > 0
    output <- data.table(distance_bin = distance_bin(distance), supported = supported)[
      !is.na(distance_bin),
      .(background_regions = .N, background_supported = sum(supported)),
      by = distance_bin
    ]
    output[, background_rate := background_supported / background_regions]
    assign(target, output, envir = background_cache)
    output
  }

  for (method in methods) {
    method_id <- method
    prediction <- copy(prediction_tables[[method]])
    total_selected <- nrow(prediction)
    prediction <- prediction[target %chin% target_union]
    if (nrow(prediction)) {
      coordinates <- regions_to_table(prediction$region)
      prediction[, `:=`(
        chrom = coordinates$chrom,
        start = coordinates$start,
        end = coordinates$end,
        midpoint = floor((coordinates$start + coordinates$end) / 2)
      )]
      prediction[, distance := {
        target_tss <- tss_table[target == .BY$target & chrom == .BY$chrom, tss]
        if (length(target_tss)) {
          vapply(midpoint, function(value) {
            min(abs(value - target_tss))
          }, numeric(1))
        } else {
          rep(Inf, .N)
        }
      }, by = .(target, chrom)]
      prediction[, distance_bin := distance_bin(distance)]
      prediction[, contact_supported := {
        ranges <- GenomicRanges::GRanges(
          chrom, IRanges::IRanges(start = start, end = end)
        )
        GenomicRanges::countOverlaps(
          ranges, support_ranges[[.BY$target]],
          ignore.strand = TRUE
        ) > 0
      }, by = target]
      matched_rows <- lapply(unique(prediction$target), function(target) {
        target_id <- target
        background <- build_background(target)
        selected <- prediction[target == target_id]
        merge(selected, background, by = "distance_bin", all.x = TRUE)
      })
      annotated <- rbindlist(matched_rows, use.names = TRUE, fill = TRUE)
      annotated <- annotated[
        is.finite(background_rate) & background_regions > 0
      ]
    } else {
      annotated <- data.table()
    }
    observed <- if (nrow(annotated)) sum(annotated$contact_supported) else 0
    expected <- if (nrow(annotated)) sum(annotated$background_rate) else 0
    variance <- if (nrow(annotated)) {
      sum(annotated$background_rate * (1 - annotated$background_rate))
    } else {
      0
    }
    z <- if (variance > 0) (observed - expected) / sqrt(variance) else NA_real_
    result_index <- result_index + 1L
    result_rows[[result_index]] <- data.frame(
      method = method,
      timepoint = timepoint,
      selected_region_gene_relations = total_selected,
      selected_with_HSPC_anchor_annotated_target = nrow(prediction),
      selected_with_distance_matched_background = nrow(annotated),
      contact_supported = observed,
      contact_support_fraction = if (nrow(annotated)) observed / nrow(annotated) else NA_real_,
      distance_matched_expected_supported = expected,
      distance_matched_expected_fraction = if (nrow(annotated)) expected / nrow(annotated) else NA_real_,
      enrichment_over_distance_matched_expectation = if (expected > 0) observed / expected else NA_real_,
      edge_level_normal_approximation_z_descriptive = z,
      edge_level_normal_approximation_p_descriptive = if (is.finite(z)) {
        stats::pnorm(z, lower.tail = FALSE)
      } else {
        NA_real_
      },
      matching_policy = "target_gene_and_prespecified_genomic_distance_bin",
      random_sampling = FALSE,
      stringsAsFactors = FALSE
    )
    if (nrow(annotated)) {
      edge_index <- edge_index + 1L
      annotated[, `:=`(method = method_id, timepoint = timepoint)]
      edge_rows[[edge_index]] <- annotated[, .(
        method, timepoint, region, target, distance, distance_bin,
        contact_supported, background_regions, background_supported,
        background_rate
      )]
    }
  }
}

edge_annotations <- rbindlist(edge_rows, use.names = TRUE, fill = TRUE)
target_gene_summary <- edge_annotations[,
  {
    observed <- sum(contact_supported)
    expected <- sum(background_rate)
    observed_fraction <- observed / .N
    expected_fraction <- expected / .N
    list(
      selected_edges = .N,
      contact_supported = observed,
      expected_supported = expected,
      observed_support_fraction = observed_fraction,
      expected_support_fraction = expected_fraction,
      support_fraction_difference = observed_fraction - expected_fraction,
      enrichment_over_distance_matched_expectation = if (expected > 0) {
        observed / expected
      } else {
        NA_real_
      }
    )
  },
  by = .(method, timepoint, target)
]

target_gene_tests <- target_gene_summary[,
  {
    differences <- support_fraction_difference[is.finite(support_fraction_difference)]
    list(
      target_gene_blocks = length(differences),
      mean_target_gene_support_difference = if (length(differences)) {
        mean(differences)
      } else {
        NA_real_
      },
      median_target_gene_support_difference = if (length(differences)) {
        stats::median(differences)
      } else {
        NA_real_
      },
      target_gene_block_wilcoxon_p = if (
        length(differences) >= 3L && any(differences != 0)
      ) {
        suppressWarnings(stats::wilcox.test(
          differences,
          mu = 0, alternative = "greater", exact = FALSE
        )$p.value)
      } else {
        NA_real_
      }
    )
  },
  by = .(method, timepoint)
]
target_gene_tests[, `:=`(
  target_gene_block_holm_p = p.adjust(
    target_gene_block_wilcoxon_p,
    method = "holm"
  ),
  target_gene_block_BH_q = p.adjust(
    target_gene_block_wilcoxon_p,
    method = "BH"
  )
)]

results <- merge(
  rbindlist(result_rows, use.names = TRUE), target_gene_tests,
  by = c("method", "timepoint"), all.x = TRUE, sort = FALSE
)
if (nrow(results) != length(methods) * length(timepoints) ||
  any(results$random_sampling)) {
  stop("HSPC region-gene contact evaluation accounting failed.")
}
write_tsv(results, file.path(stage_root, "HSPC_contact_summary.tsv"))
write_tsv(edge_annotations, file.path(stage_root, "HSPC_contact_edge_annotations.tsv.gz"))
write_tsv(
  target_gene_summary,
  file.path(stage_root, "HSPC_contact_target_gene_summary.tsv.gz")
)
write_tsv(support_table, file.path(stage_root, "HSPC_loop_anchor_gene_map.tsv.gz"))

if (!file.rename(stage_root, output_root)) {
  stop("Unable to atomically publish HSPC contact evaluation.")
}
cat("region_gene_HSPC_contact_evaluation\tPASS\t", output_root, "\n", sep = "")

}

if (sys.nframe() == 0L) evaluate_region_gene_contact()
