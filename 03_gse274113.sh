#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source code/functions/runtime.sh
run_task code/gse274113/build_control_object.R
run_task code/gse274113/prepare_control_candidates.R
run_task code/gse274113/fit_control_timepoints.R
run_task code/gse274113/fit_replicate_stability.R
run_task code/gse274113/summarize_leave_one_replicate.R
run_task code/gse274113/summarize_pairwise_fold_stability.R
run_task code/gse274113/summarize_domain_decomposition.R
run_task code/gse274113/summarize_domain_stability.R
run_task code/gse274113/build_heldout_truth.R
run_task code/gse274113/run_pando.R
run_task code/gse274113/prepare_scenicplus_inputs.py
run_task code/gse274113/run_pycistopic.py
run_task code/gse274113/run_scenicplus.py
run_task code/gse274113/export_endpoints.R
run_task code/gse274113/summarize_network_catalog.R
run_task code/gse274113/evaluate_tf_gene_comparators.R
run_task code/gse274113/summarize_method_evidence.R
