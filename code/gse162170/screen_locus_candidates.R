#!/usr/bin/env Rscript


suppressPackageStartupMessages(library(data.table))

# ---- Parameters ----
triplet_path <- normalizePath("results/gse162170/network/triplets.tsv", mustWork = TRUE)
membership_path <- normalizePath("results/gse162170/cell_programs/target_membership.tsv", mustWork = TRUE)
output_dir <- "results/gse162170/loci/candidate_screen"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

triplets <- fread(triplet_path, showProgress = FALSE)
membership <- fread(membership_path)
summary <- triplets[, .(
  chains = .N,
  regions = uniqueN(region),
  regulators = uniqueN(regulator),
  positive = sum(chain_direction > 0),
  negative = sum(chain_direction < 0),
  median_chain_delta_bic = median(chain_delta_bic),
  q75_chain_delta_bic = quantile(chain_delta_bic, 0.75),
  max_chain_delta_bic = max(chain_delta_bic),
  backbone_selected = sum(backbone_selected),
  backbone_sign_concordant = sum(backbone_sign_concordant)
), by = target]
summary <- merge(summary, membership, by = "target", all.x = TRUE, sort = FALSE)
summary[, sign_balance := pmin(positive, negative) / pmax(positive, negative)]
summary[, eligible :=
  rna_matches_chromatin_anchor %in% TRUE &
    regions >= 5L & regions <= 15L &
    chains >= 20L & chains <= 60L &
    positive >= 5L & negative >= 5L &
    median_chain_delta_bic >= 3 & target_margin >= 0.1]
setorder(
  summary, -eligible, -sign_balance, -median_chain_delta_bic,
  -target_margin, -chains, target
)
summary[, evidence_rank := NA_integer_]
summary[eligible == TRUE, evidence_rank := seq_len(.N)]

fwrite(summary, file.path(output_dir, "all_target_locus_screen.tsv"), sep = "\t")
fwrite(summary[eligible == TRUE], file.path(output_dir, "eligible_locus_shortlist.tsv"), sep = "\t")
