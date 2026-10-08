#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source code/functions/runtime.sh
run_task code/gse162170/fit_network.R
run_task code/gse162170/extract_cell_programs.R
run_task code/gse162170/summarize_cell_programs.R
run_task code/gse162170/analyze_key_gene_gsea.R
run_task code/gse162170/run_gsea.R
run_task code/gse162170/extract_locus_evidence.R
run_task code/gse162170/materialize_locus_evidence.R
run_task code/gse162170/screen_locus_candidates.R
run_task code/gse162170/rank_loci_by_cell_type.R
