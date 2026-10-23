#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source code/functions/runtime.sh
run_task code/gse162170/prepare_candidates.R
run_task code/gse162170/fit_network.R
run_task code/gse162170/extract_cell_programs.R
run_task code/gse162170/summarize_aucell_thresholds.R
run_task code/gse162170/summarize_cell_programs.R
run_task code/gse162170/analyze_key_gene_gsea.R
run_task code/gse162170/run_gsea.R
run_task code/gse162170/extract_locus_evidence.R
run_task code/gse162170/materialize_locus_evidence.R
run_task code/gse162170/screen_locus_candidates.R
run_task code/gse162170/rank_loci_by_cell_type.R
run_task code/gse162170/evaluate_locus_ordering.py
run_task code/methods/fit_cell_budget_networks.R run_cortex_sensitivity
run_task code/gse162170/summarize_within_class_sensitivity.R
