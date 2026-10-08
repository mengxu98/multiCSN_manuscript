#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Assemble the Source Data workbook from the formal result tables.

Every sheet is built from formal results (including shared BEELINE case views), so that the Source Data file
shipped with the manuscript stays in step with the analysis outputs. The index
sheet describes each sheet by figure panel and content only; internal file
locations and intermediate artifacts are deliberately not recorded there.

Supplementary Tables S1 and S2 are typeset in the supplementary PDF and are
deliberately not repeated here.
"""

import os
import re
import sys
import subprocess
import tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "functions"))
from runtime import use_python, task_environment, r_function_command
if __name__ == "__main__":
    use_python("general")

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

PROJECT_ROOT = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)
OUT_DIR = os.path.join(PROJECT_ROOT, "source_data")
CSV_DIR = os.path.join(OUT_DIR, "files")
SOURCE_XLSX = os.path.join(OUT_DIR, "Source Data.xlsx")
FONT_NAME = "Times New Roman"
HEADER_FILL = "D3D3D3"  # light gray
MAX_SHEET_NAME = 31

# Column headings keep these tokens in upper case.
_PRESERVE_UPPER = {
    "TF", "TFS", "CRE", "CRES", "RNA", "ATAC", "DNA", "GO", "BH", "FDR", "NES",
    "AP", "AUPRC", "AUROC", "EPR", "HSPC", "GSEA", "AUC", "AUCell", "UMAP",
    "PC", "PCA", "BIC", "RSS", "SNR", "ID", "IDS", "Jaccard", "Wilcoxon",
    "Holm", "KCNC2", "ZNF804A", "MK/E", "AUCELL", "SE",
}
_PRESERVE_MIXED = {
    "AUCell": "AUCell", "AUPRC": "AUPRC", "AUROC": "AUROC", "NES": "NES",
    "EPR": "EPR", "HSPC": "HSPC", "GSEA": "GSEA", "Jaccard": "Jaccard",
    "Wilcoxon": "Wilcoxon", "Holm": "Holm", "Benjamini": "Benjamini",
    "Hochberg": "Hochberg", "TF": "TF", "TFs": "TFs", "CRE": "CRE",
    "CRES": "CREs", "RNA": "RNA", "ATAC": "ATAC", "DNA": "DNA", "GO": "GO",
    "BH": "BH", "FDR": "FDR", "AP": "AP", "AUC": "AUC", "BIC": "BIC",
    "UMAP": "UMAP",
}


def log(message, level="info"):
    prefix = {"info": "[info]", "warn": "[warn]", "fail": "[fail]", "success": "[ok]"}.get(level, f"[{level}]")
    print(f"{prefix} {message}", flush=True)


def _format_column_name(name):
    """Turn a result-table column name into a readable heading."""
    key = str(name).strip()
    if key in _PRESERVE_MIXED:
        return _PRESERVE_MIXED[key]
    words = re.split(r"[_\s]+", key)
    out = []
    for i, w in enumerate(words):
        if not w:
            continue
        if w in _PRESERVE_MIXED:
            out.append(_PRESERVE_MIXED[w])
        elif w.upper() in _PRESERVE_UPPER:
            out.append(w.upper())
        elif w.lower() in ("of", "and", "per", "vs", "by", "in", "to") and i:
            out.append(w.lower())
        else:
            out.append(w[:1].upper() + w[1:])
    return " ".join(out) if out else key


def _format_columns(df):
    df = df.copy()
    df.columns = [_format_column_name(c) for c in df.columns]
    return df


def _sanitize_sheet_name(name):
    clean = re.sub(r"[\[\]:*?/\\]", "-", name).strip()
    if len(clean) > MAX_SHEET_NAME:
        clean = clean[:MAX_SHEET_NAME]
    return clean


def read_table(relative_path):
    """Read a result table (tsv/csv, optionally gzip-compressed)."""
    path = os.path.join(PROJECT_ROOT, relative_path)
    if not os.path.exists(path):
        raise FileNotFoundError(relative_path)
    if path.endswith(".gz"):
        engine_path = path
        lowered = path[:-3].lower()
    else:
        engine_path = path
        lowered = path.lower()
    sep = "," if lowered.endswith(".csv") else "\t"
    return pd.read_csv(engine_path, sep=sep, compression="infer", low_memory=False)


# --------------------------------------------------------------------------
# Source Data sheets: one entry per formal result table.
#   figure  : figure panel(s) the table supports
#   content : what the table contains, in article terms
# --------------------------------------------------------------------------
SOURCE_SHEETS = [
    # Figure 2 - case-matched BEELINE evaluation
    dict(sheet="Fig2a_case_metrics", figure="Figure 2a",
         content="Case-level AUPRC, EarlyPrecision and EPR for 500 BEELINE cases and 11 methods",
         path="results/beeline/metrics/case_metrics.tsv"),
    dict(sheet="Fig2b_rank_summary", figure="Figure 2b",
         content="Mean, median and quartile within-case ranks per method and metric across the 500 cases",
         path="results/beeline/metrics/rank_summary.tsv"),
    dict(sheet="Fig2b_case_ranks", figure="Figure 2b",
         content="Within-case rank of every method and metric in every BEELINE case",
         path="results/beeline/metrics/case_ranks.tsv"),
    dict(sheet="Fig2c_runtime_summary", figure="Figure 2c",
         content="Single-worker inference-and-export runtime totals, medians, maxima and ratios to multiCSN (500 cases per method; input reading and matrix preparation excluded)",
         path="results/beeline/runtime/runtime_summary.tsv"),
    dict(sheet="Fig2c_runtime_all_methods", figure="Figure 2c",
         content="Per-case inference-and-export runtimes for 11 methods; input reading and matrix preparation are excluded and recorded separately as LoadSeconds in runtime/case_runtime.tsv",
         path="results/beeline/runtime/runtime_all_methods.tsv"),
    dict(sheet="Fig2d_HSC_truth", figure="Figure 2d",
         content="Reference network for the HSC 2000-3-30 BEELINE case",
         case="HSC-2000-3-30", view="truth"),
    dict(sheet="Fig2d_HSC_methods", figure="Figure 2d",
         content="Per-method edge counts, precision, recall, F1 and sign concordance for the HSC 2000-3-30 case",
         case="HSC-2000-3-30", view="summary"),
    dict(sheet="Fig2d_HSC_edges", figure="Figure 2d",
         content="Displayed edges of every method for the HSC 2000-3-30 case at the requested edge budget",
         case="HSC-2000-3-30", view="edges"),

    dict(sheet="Resource_code_provenance", figure="Supplementary Table S2",
         content="Code pins, archived source hashes and verified motif resource provenance",
         path="source_data/supplementary_tables/resource_and_code_provenance.tsv"),

    # Figure 3 - hematopoietic Perturb-multiome
    dict(sheet="Fig3a_response_metrics", figure="Figure 3a",
         content="Equal-TF mean paired differences, TF-bootstrap intervals and adjusted P values",
         path="source_data/files/fig3a_multimetric_effects_source.tsv"),
    dict(sheet="Fig3a_TF_paired_effects", figure="Figure 3a",
         content="TF-level mean paired differences and finite timepoint counts underlying equal-TF effects and intervals",
         path="source_data/files/fig3a_tf_cluster_effects_source.tsv"),
    dict(sheet="Fig3b_TF_region", figure="Figure 3b",
         content="UniBind TF-region precision and recall in the four-method common non-empty strata",
         path="source_data/files/fig3b_tf_region_source.tsv"),
    dict(sheet="Fig3c_region_gene", figure="Figure 3c",
         content="Target-blocked region-gene contact enrichment over distance-matched expectation",
         path="source_data/files/fig3c_region_gene_source.tsv"),
    dict(sheet="Fig3d_direct_enrichment", figure="Figure 3d",
         content="Direct RNA- and ATAC-target response enrichment per method, regulator and time point",
         path="source_data/files/fig3d_direct_enrichment_source.tsv.gz"),
    dict(sheet="Fig3ef_RUNX1_paths", figure="Figure 3e,f",
         content="RUNX1 day-9 direction-consistent direct and second-order paths",
         path="results/gse274113/perturbation_neighborhood/hero_direction_consistent_paths.tsv"),
    dict(sheet="Fig3ef_chain_annotations", figure="Figure 3e,f",
         content="Held-out RNA and ATAC responses annotated on every selected day-9 chain",
         path="results/gse274113/chain_response_annotations/selected_chain_response_annotations.tsv.gz"),
    dict(sheet="Fig3g_response_order", figure="Figure 3g",
         content="Direct and second-order RNA-response enrichment per perturbed regulator",
         path="source_data/files/fig3g_order_response_source.tsv"),
    dict(sheet="Fig3h_concordance", figure="Figure 3h",
         content="Signed response concordance per method, regulator and time point",
         path="source_data/files/fig3h_response_concordance_source.tsv"),
    dict(sheet="Fig3h_summary", figure="Figure 3h",
         content="Median and interquartile signed response concordance per method and time point",
         path="source_data/files/fig3h_response_concordance_summary.tsv"),
    dict(sheet="Fig3i_layer_attrition", figure="Figure 3i",
         content="Replicate-deletion retention of complete-network edges by eligibility layer and fold",
         path="results/gse274113/replicate_deletion_stability/multilayer_decomposition/layer_attrition.tsv"),
    dict(sheet="Fig3j_chain_response", figure="Figure 3j",
         content="Joint ATAC-RNA response support of selected chains per method and TF cluster",
         path="source_data/files/fig3j_chain_response_source.tsv"),

    # Figure 4 - developing human fetal cortex
    dict(sheet="Fig4a_topology", figure="Figure 4a",
         content="Candidate-pair and selected-edge counts for the three component networks",
         path="results/gse162170/network/figure_source/network_topology.tsv"),
    dict(sheet="Fig4b_target_funnel", figure="Figure 4b",
         content="Target funnel from annotated genes to selected sparse support",
         path="results/gse162170/network/figure_source/Fig4b_target_funnel.tsv"),
    dict(sheet="Fig4c_chain_summary", figure="Figure 4c",
         content="Exact-region chain counts and direct TF-gene overlap",
         path="results/gse162170/network/figure_source/chain_summary.tsv"),
    dict(sheet="Fig4d_terminal_states", figure="Figure 4d",
         content="Mutually exclusive terminal target classes",
         path="results/gse162170/network/figure_source/terminal_states.tsv"),
    dict(sheet="Fig4e_aucell_by_class", figure="Figure 4e",
         content="RNA AUCell of the five chromatin-anchored target programs by cell class",
         path="results/gse162170/cell_programs/aucell_by_class.tsv"),
    dict(sheet="Fig4e_aucell_by_cell", figure="Figure 4e",
         content="Per-cell RNA AUCell scores, target-program identity and cell class",
         path="results/gse162170/cell_programs/aucell_by_cell.tsv.gz",
         expected_rows=44905),
    dict(sheet="Fig4e_threshold_summary", figure="Figure 4e",
         content="Fixed program-specific AUC thresholds, above-threshold counts, cell totals and fractions by cell class",
         path="results/gse162170/cell_programs/aucell_threshold_summary.tsv",
         expected_rows=25),
    dict(sheet="Fig4f_module_activity", figure="Figure 4f",
         content="Aligned RNA target-AUCell and selected-region ATAC activity for the 25 regulator-anchor-sign modules",
         path="results/gse162170/cell_programs/module_raw_and_scaled_activity.tsv"),
    dict(sheet="Fig4f_S8_module_summary", figure="Figure 4f; Figure S8",
         content="Membership and size summary of the 25 regulator-anchor-sign modules",
         path="results/gse162170/cell_programs/module_summary.tsv"),
    dict(sheet="Fig4g_go_terms", figure="Figure 4g",
         content="GO biological-process enrichment of the five target-program families",
         path="results/gse162170/cell_programs/go_display_terms.tsv"),
    dict(sheet="Fig4h_gsea_terms", figure="Figure 4h",
         content="Displayed preranked GSEA terms over the selected-target RNA rank",
         path="results/gse162170/key_gene_gsea/gsea_selected_terms.tsv"),
    dict(sheet="Fig4h_gsea_genes", figure="Figure 4h",
         content="Leading-edge and rank statistics of the displayed GSEA terms",
         path="results/gse162170/key_gene_gsea/gsea_key_genes.tsv"),
    dict(sheet="Fig4h_rank_curve", figure="Figure 4h",
         content="Complete 11,017-target IPC/RG RNA ranking, exact plotted running enrichment score, pathway membership and leading-edge membership",
         path="results/gse162170/key_gene_gsea/gsea_rank_curve.tsv.gz",
         expected_rows=11017),
    dict(sheet="Fig4ij_5_ZNF804A_paths", figure="Figure 4i,j; Figure S10a,b",
         content="The 32 ZNF804A exact-region paths with regulators, regions, signs and deletion evidence",
         path="results/gse162170/loci/candidates/ZNF804A_paths.tsv",
         identity=("target", "ZNF804A"), expected_rows=32),
    dict(sheet="Fig4i_ZNF804A_annotation", figure="Figure 4i",
         content="ZNF804A gene annotation used with the path coordinates to draw the locus track",
         path="results/gse162170/loci/candidates/ZNF804A_gene_annotation.tsv",
         identity=("gene_name", "ZNF804A")),
    dict(sheet="FigS10a_ZNF804A_expression", figure="Figure S10a",
         content="Per-sample ZNF804A RNA expression by cell class",
         path="results/gse162170/loci/candidates/ZNF804A_expression_samples.tsv",
         identity=("target", "ZNF804A")),
    dict(sheet="Fig4i_ZNF804A_accessibility", figure="Figure 4i",
         content="Region accessibility by cell class for the ZNF804A locus track",
         path="results/gse162170/loci/candidates/ZNF804A_accessibility.tsv"),
    dict(sheet="FigS10_KCNC2_paths", figure="Figure S10c-e",
         content="The 33 KCNC2 exact-region paths with regulators, regions, signs and deletion evidence",
         path="results/gse162170/loci/candidates/KCNC2_paths.tsv",
         identity=("target", "KCNC2"), expected_rows=33),
    dict(sheet="FigS10c_KCNC2_annotation", figure="Figure S10c",
         content="KCNC2 gene annotation used with the path coordinates to draw the locus track",
         path="results/gse162170/loci/candidates/KCNC2_gene_annotation.tsv",
         identity=("gene_name", "KCNC2")),
    dict(sheet="FigS10d_KCNC2_expression", figure="Figure S10d",
         content="Per-sample KCNC2 RNA expression by cell class",
         path="results/gse162170/loci/candidates/KCNC2_expression_samples.tsv",
         identity=("target", "KCNC2")),
    dict(sheet="FigS10c_KCNC2_accessibility", figure="Figure S10c",
         content="Region accessibility by cell class for the KCNC2 locus track",
         path="results/gse162170/loci/candidates/KCNC2_accessibility.tsv"),

    # Figure 5 - human erythroid maturation
    dict(sheet="Fig5a_cells_pseudotime", figure="Figure 5a",
         content="Gene-expression UMAP coordinates and pseudotime of the 10,862 erythroid cells with finite pseudotime",
         path="results/gse194122/trajectory/erythroid_coordinates_and_pseudotime.tsv.gz"),
    dict(sheet="Fig5b_g_focal_regulons", figure="Figure 5b,g",
         content="The 582 regulons with focal-regulon selection statistics and HSC-to-Normoblast activity change",
         path="results/gse194122/state_regulomes/hero_regulons.tsv"),
    dict(sheet="Fig5b_regulon_network", figure="Figure 5b",
         content="Global map of the six focal positive regulons with donor-stable chain support",
         path="results/gse194122/state_regulomes/hero_regulon_network.tsv.gz"),
    dict(sheet="Fig5c_state_activity", figure="Figure 5c",
         content="State-aggregated gene-module expression and CRE-linked network support per regulon and state",
         path="results/gse194122/state_regulomes/module_state_summary.tsv.gz"),
    dict(sheet="Fig5d_feature_programs", figure="Figure 5d",
         content="Temporal cluster assignment of the 446 focal TF and target genes",
         path="source_data/files/dynamic_feature_programs.tsv"),
    dict(sheet="Fig5d_go_enrichment", figure="Figure 5d",
         content="GO biological-process enrichment per temporal cluster of the fitted trajectories",
         path="source_data/files/dynamic_go_enrichment.tsv.gz"),
    dict(sheet="Fig5e_h_model_variance", figure="Figure 5e,h",
         content="Per-regulon activity coupling and partial state and donor R-squared values",
         path="results/gse194122/state_regulomes/module_model_variance.tsv.gz"),
    dict(sheet="Fig5f_transition_topology", figure="Figure 5f",
         content="Gained, lost, sign-stable and sign-flipped support across adjacent erythroid states",
         path="results/gse194122/state_transitions/transition_summary/GSE5d_adjacent_transition_topology.tsv"),
    dict(sheet="Fig5ij_rewiring_GO", figure="Figure 5i,j",
         content="BH-significant GO biological-process terms of donor-stable target modules gained and lost across transitions",
         path="results/gse194122/donor_stable_rewiring/stable_rewiring_GO_BP.tsv.gz",
         filter=("p_adjust", "<", 0.05)),

    # Supplementary Figures
    dict(sheet="FigS1_BEELINE_other_metrics", figure="Figure S1",
         content="AUROC, accuracy, precision, recall and F1 per case and method across the 500 BEELINE cases",
         path="results/beeline/metrics/combined__other_metrics.csv"),
    dict(sheet="FigS2_curated_dropout", figure="Figure S2",
         content="AUPRC, EarlyPrecision and EPR for the 200 curated-model cases across dropout conditions",
         path="results/beeline/metrics/curated__primary_metrics.csv"),
    dict(sheet="FigS3_synthetic_cells", figure="Figure S3",
         content="AUPRC, EarlyPrecision and EPR for the 300 synthetic-network cases across cell-number settings",
         path="results/beeline/metrics/synthetic__primary_metrics.csv"),
    dict(sheet="FigS4_scaling_runs", figure="Figure S4",
         content="Runtime and sampled process-tree RSS sum for 60 fits across cell counts, total gene counts and requested cores",
         path="results/beeline/core_scaling/core_scaling_runs.csv",
         columns=["Case", "Scenario", "AxisValue", "AxisLabel", "Cores", "InputCells", "InputGenes",
                  "RuntimeSeconds", "UserSeconds", "SystemSeconds", "OutputEdges",
                  "MemoryTreePeakMb", "MemoryTreeBeforeMb", "MemoryTreeProcessesAtPeak", "MemorySamples",
                  "MemorySamplingIntervalSeconds", "MemoryMaxObservedIntervalSeconds",
                  "TimingScope"]),
    dict(sheet="FigS4_core_comparison", figure="Figure S4",
         content="One-core versus ten-core inference runtime and peak RSS sum (MiB) per scaling configuration with selected-edge counts",
         path="results/beeline/core_scaling/one_to_ten_core_comparison.csv"),
    dict(sheet="FigS4_scaling_summary", figure="Figure S4",
         content="Summary speedup and memory ratios of the scaling experiment",
         path="results/beeline/core_scaling/scaling_summary_statistics.csv"),
    dict(sheet="FigS4_edge_invariance", figure="Figure S4",
         content="Selected-edge counts across requested core settings; this table does not establish full-network equality",
         path="results/beeline/core_scaling/edge_count_invariance.csv",
         columns=["Case", "CoreSettings", "EdgeCounts", "OutputEdges"]),
    dict(sheet="FigS5_VSC_truth", figure="Figure S5",
         content="Reference network for the VSC 2000-9 BEELINE case",
         case="VSC-2000-9", view="truth"),
    dict(sheet="FigS5_VSC_methods", figure="Figure S5",
         content="Per-method recovery and sign concordance for the VSC 2000-9 case",
         case="VSC-2000-9", view="summary"),
    dict(sheet="FigS5_VSC_edges", figure="Figure S5",
         content="Displayed edges of every method for the VSC 2000-9 case",
         case="VSC-2000-9", view="edges"),
    dict(sheet="FigS6_scMultiSim_seeds", figure="Figure S6a-c",
         content="Per-seed TF-gene and region-gene recovery metrics across the 10 scMultiSim simulations",
         path="results/scmultisim/multiseed/endpoint_metrics_all_seeds.tsv"),
    dict(sheet="FigS6_scMultiSim_diagnostics", figure="Figure S6a-c",
         content="Per-seed ranking diagnostics and prevalence baselines",
         path="results/scmultisim/multiseed_analysis/per_seed_diagnostics.tsv"),
    dict(sheet="FigS6_truth_sensitivity", figure="Figure S6c",
         content="Paired static-input-GRN versus dynamic-GRN-union recovery per seed",
         path="results/scmultisim/multiseed_analysis/static_dynamic_truth_sensitivity.tsv"),
    dict(sheet="FigS6_ablation_seeds", figure="Figure S6d-f",
         content="Per-seed TF-gene recovery for the forward-only and full add-delete-swap search variants",
         path="results/scmultisim/component_ablation/endpoint_metrics_all_seeds.tsv"),
    dict(sheet="FigS7_runtime", figure="Figure S7a (runtime)",
         content="Inference-only runtime of the four multi-omic methods on GSE274113, with the measured core count and scope of each configuration",
         path="results/gse274113/method_comparison/multiomic_runtime_comparison.tsv"),
    dict(sheet="FigS7_fold_stability", figure="Figure S7c",
         content="Per-fold replicate-deletion stability metrics",
         path="results/gse274113/replicate_deletion_stability/multilayer_decomposition/all_fold_metrics.tsv"),
    dict(sheet="FigS7_timepoint_summary", figure="Figure S7 (stability summary)",
         content="Per-time-point medians and quartiles of replicate-deletion stability metrics",
         path="results/gse274113/replicate_deletion_stability/summary/timepoint_stability_summary.tsv"),
    dict(sheet="FigS7_pairwise_stability", figure="Figure S7c",
         content="Edge-set Jaccard similarity for fold pairs and fold-versus-complete comparisons",
         path="results/gse274113/replicate_deletion_stability/pairwise_stability/pairwise_fold_metrics.tsv"),
    dict(sheet="FigS7_rank_deciles", figure="Figure S7e",
         content="Edge retention by within-network rank decile",
         path="results/gse274113/replicate_deletion_stability/multilayer_decomposition/rank_decile_decomposition.tsv"),
    dict(sheet="FigS7_top100", figure="Figure S7d",
         content="Top-100 edge rank decomposition of replicate-deletion retention",
         path="results/gse274113/replicate_deletion_stability/summary/top100_decomposition.tsv"),
    dict(sheet="FigS7e_endpoint_counts", figure="Figure S7b",
         content="Selected TF-gene, TF-region, region-gene and exact-chain counts at each time point",
         path="results/gse274113/factorized_chains/timepoint_summary.tsv"),
    dict(sheet="FigS8_module_membership", figure="Figure S8",
         content="Target membership of the 25 regulator-anchor-sign modules",
         path="results/gse162170/cell_programs/module_membership.tsv.gz"),
    dict(sheet="FigS11_regulator_turnover", figure="Figure S11",
         content="Tie-inclusive TF-gene and factorized-chain regulator turnover across adjacent erythroid states",
         path="results/gse194122/state_transitions/transition_summary/GSE5e_regulator_transition_counts_display.tsv"),
    dict(sheet="FigS12_cell_budget_a", figure="Figure S12a",
         content="Signed Jaccard similarity between each fixed-cell-budget fit and the all-cell fit",
         path="source_data/files/figs12_panel_a_state_overlap.tsv"),
    dict(sheet="FigS12_cell_budget_b", figure="Figure S12b",
         content="Adjacent-state topology-composition similarity under fixed-cell-budget refitting",
         path="source_data/files/figs12_panel_b_transition_similarity.tsv"),
    dict(sheet="FigS12_cell_budget_c", figure="Figure S12c",
         content="Focal-regulon target support in the fixed-cell-budget fits versus the all-cell fits",
         path="source_data/files/figs12_panel_c_focal_module_support.tsv"),
    dict(sheet="FigS12_cell_budget_d", figure="Figure S12d",
         content="Selected exact-region chain counts divided by fixed consensus-module target counts",
         path="source_data/files/figs12_panel_d_focal_chain_support.tsv"),
    dict(sheet="TableS3_BEELINE_means", figure="Table S3",
         content="Case-weighted AUPRC means and mean within-case ranks across the same 500 cases",
         path="results/beeline/evaluation/beeline_method_summary.tsv"),
    dict(sheet="BEELINE_case_manifest", figure="Table S3",
         content="Model, cells, dropout, replicate and gene/truth-edge counts from the reference networks",
         path="results/beeline/evaluation/beeline_case_manifest.tsv"),
    dict(sheet="TableS4_budget_summary", figure="Table S4",
         content="Focal-module peak and expression-support correlation sensitivity; RNA summaries held fixed",
         path="results/gse194122/cell_budget_sensitivity/focal_budget_summary.tsv"),
    dict(sheet="TableS4_budget_draws", figure="Table S4",
         content="Per-draw peak sets and correlations for target and exact-chain support",
         path="results/gse194122/cell_budget_sensitivity/focal_budget_draws.tsv"),
    dict(sheet="FigS9a_b_locus_screen", figure="Figure S9a,b",
         content="Complete ordered fetal-cortex locus screen with all screening statistics and within-class ranks",
         path="results/gse162170/loci/by_class/eligible_loci_by_class.tsv"),
]



def _apply_filter(df, spec):
    column, op, value = spec
    if column not in df.columns:
        raise KeyError(f"filter column {column!r} not in table")
    if op == "<":
        return df[df[column] < value]
    if op == "<=":
        return df[df[column] <= value]
    if op == ">":
        return df[df[column] > value]
    if op == ">=":
        return df[df[column] >= value]
    if op == "==":
        return df[df[column] == value]
    raise ValueError(f"unsupported filter operator {op!r}")


_case_views = {}


def load_source_sheet(entry):
    if "case" in entry:
        case = entry["case"]
        if case not in _case_views:
            with tempfile.TemporaryDirectory(prefix="beeline-case-") as directory:
                subprocess.run(r_function_command(
                    "code/functions/beeline_case_views.R", "write_beeline_case", role="beeline",
                    repo=PROJECT_ROOT, case=case, output=directory), env=task_environment(), check=True)
                _case_views[case] = {
                    name: pd.read_csv(os.path.join(directory, name + ".tsv"), sep="\t")
                    for name in ("truth", "edges", "summary")
                }
        df = _case_views[case][entry["view"]].copy()
    else:
        df = read_table(entry["path"])
    if entry.get("identity"):
        column, expected = entry["identity"]
        if column not in df or df.empty or df[column].isna().any() or set(df[column]) != {expected}:
            raise ValueError(f"{entry['sheet']}: expected {column}={expected} in {entry['path']}")
    if "expected_rows" in entry and len(df) != entry["expected_rows"]:
        raise ValueError(f"{entry['sheet']}: expected {entry['expected_rows']} rows, got {len(df)}")
    if entry.get("filter"):
        df = _apply_filter(df, entry["filter"])
    if entry.get("columns"):
        df = df[entry["columns"]]
    return _format_columns(df)


def write_workbook(xlsx_path, index, sheets, csv_prefix):
    """Write the index sheet, one sheet per table, and per-sheet CSV copies."""
    used = set()
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        index.to_excel(writer, sheet_name="Table information", index=False)
        for name, df in sheets:
            sheet_name = _sanitize_sheet_name(name)
            base = sheet_name
            suffix = 1
            while sheet_name in used:
                suffix += 1
                sheet_name = _sanitize_sheet_name(f"{base}_{suffix}")
            used.add(sheet_name)
            df.to_excel(writer, sheet_name=sheet_name, index=False)
            df.to_csv(os.path.join(CSV_DIR, f"{csv_prefix}{sheet_name}.csv"), index=False)
            log(f"{sheet_name}: {len(df):,} rows")

    wb = load_workbook(xlsx_path)
    header_fill = PatternFill(fill_type="solid", fgColor=HEADER_FILL)
    header_font = Font(name=FONT_NAME, bold=True)
    body_font = Font(name=FONT_NAME)
    for ws in wb.worksheets:
        for col in range(1, ws.max_column + 1):
            cell = ws.cell(row=1, column=col)
            cell.fill = header_fill
            cell.font = header_font
        for row in ws.iter_rows(min_row=2, max_row=ws.max_row, max_col=ws.max_column):
            for cell in row:
                cell.font = body_font
        for col in range(1, ws.max_column + 1):
            width = 0
            for row in range(1, min(ws.max_row + 1, 102)):
                value = ws.cell(row=row, column=col).value
                if value is not None:
                    width = max(width, len(str(value)))
            ws.column_dimensions[get_column_letter(col)].width = min(max(width + 1, 10), 50)
        if ws.title in {"Resource_code_provenance", "Fig3a_TF_paired_effects"}:
            ws.freeze_panes = "A2"
            ws.row_dimensions[1].height = 32
            for cell in ws[1]:
                cell.alignment = Alignment(wrap_text=True, vertical="center")
            for row in ws.iter_rows(min_row=2):
                for cell in row:
                    cell.alignment = Alignment(wrap_text=True, vertical="top")
                ws.row_dimensions[row[0].row].height = 110 if ws.title == "Resource_code_provenance" else 30
    wb.save(xlsx_path)


def main():
    _case_views.clear()
    os.makedirs(CSV_DIR, exist_ok=True)
    log(f"Building Source Data from {len(SOURCE_SHEETS)} result tables")

    index_rows = []
    sheets = []
    missing = []
    for entry in SOURCE_SHEETS:
        try:
            df = load_source_sheet(entry)
        except FileNotFoundError:
            missing.append(entry.get("path", entry.get("case")))
            log(f"{entry['sheet']}: source table not found", "fail")
            continue
        sheets.append((entry["sheet"], df))
        index_rows.append({
            "Figure": entry["figure"],
            "Sheet": _sanitize_sheet_name(entry["sheet"]),
            "Content": entry["content"],
            "Rows": len(df),
            "Columns": df.shape[1],
        })

    if missing:
        log(f"{len(missing)} source tables are missing; fix the list before shipping", "fail")
        for path in missing:
            log(f"  missing: {path}", "fail")
        return 1

    source_index = pd.DataFrame(index_rows).sort_values(
        "Figure", key=lambda s: s.map(_figure_sort_key)
    ).reset_index(drop=True)
    write_workbook(
        SOURCE_XLSX,
        source_index,
        sheets,
        csv_prefix="",
    )
    log(f"Source Data written: {SOURCE_XLSX}", "info")

    return 0


def _figure_sort_key(label):
    """Sort figure labels as Figure 2 < Figure 3 < ... < Figure S1 < Figure S10."""
    text = str(label)
    main = re.match(r"Figure (\d+)", text)
    if main:
        return (0, int(main.group(1)), 0, text)
    supp = re.search(r"Figure S(\d+)", text)
    if supp:
        return (1, int(supp.group(1)), 0, text)
    return (2, 0, 0, text)


if __name__ == "__main__":
    sys.exit(main())
