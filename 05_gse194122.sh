#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source code/functions/runtime.sh
run_task code/requirements/prepare_go_resources.R
run_task code/gse194122/prepare_candidates.R
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
run_task code/gse194122/summarize_focal_cell_budget.py
run_task code/gse194122/export_go_annotation.R
run_task code/methods/fit_cell_budget_networks.R run_erythroid_sensitivity
run_task code/gse194122/summarize_donor_exclusion_sensitivity.R
run_task code/gse194122/validate_application_cases.R
run_task code/gse194122/evaluate_chain_modules.R
run_task code/gse194122/evaluate_regulatory_programs.py
run_task code/gse194122/summarize_link_recurrence.py
run_task code/gse194122/prepare_functional_programs.py
