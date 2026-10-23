#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source code/functions/runtime.sh
run_task code/gse274113/inventory_release.R
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
run_task code/gse274113/run_linger.py
run_task code/gse274113/export_endpoints.R
run_task code/gse274113/summarize_network_catalog.R
run_task code/gse274113/export_factorized_chains.R
run_task code/gse274113/evaluate_formal_chain_responses.R
run_task code/gse274113/evaluate_tf_gene_comparators.R
run_task code/gse274113/summarize_method_evidence.R
run_task code/gse274113/evaluate_linger_tf_gene.R
run_task code/gse274113/compare_linger_vs_multicsn.R
run_task code/gse274113/evaluate_tf_region_binding.R
run_task code/gse274113/evaluate_region_gene_contact.R
run_task code/gse274113/evaluate_linger_multilayer.R run_linger_binding
run_task code/gse274113/evaluate_linger_multilayer.R run_linger_contact
run_task code/gse274113/evaluate_linger_multilayer.R run_linger_chain
run_task code/gse274113/analyze_perturbation_neighborhood.R
run_task code/gse274113/evaluate_linger_perturbation.R
run_task code/gse274113/aggregate_linger_comparison.R
run_task code/gse274113/analyze_weight_stability.R
run_task code/gse274113/analyze_eligibility_conditioned_stability.R
run_task code/gse274113/evaluate_layer_ablation.R
run_task code/gse274113/summarize_runtime.py
