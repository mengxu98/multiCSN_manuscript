#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source code/functions/runtime.sh
run_task code/gse194122/fit_state_networks.R
run_task code/gse194122/summarize_state_transitions.R
run_task code/gse194122/extract_trajectory.R
run_task code/gse194122/recompute_erythroid_gex_umap.R
run_task code/gse194122/extract_donor_state_expression.R
run_task code/gse194122/prepare_donor_deletion_manifests.R
run_task code/gse194122/run_network_batch.R run_donor_deletion
run_task code/gse194122/summarize_donor_stable_rewiring.R
run_task code/gse194122/analyze_state_regulomes.R
run_task code/gse194122/prepare_cell_budget_manifests.R
run_task code/gse194122/run_network_batch.R run_cell_budget
run_task code/gse194122/summarize_cell_budget_sensitivity.R
