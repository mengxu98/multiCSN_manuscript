from pathlib import Path
import sys
import tempfile
import warnings

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "functions"))
from runtime import REPO, use_python
if __name__ == "__main__":
    use_python("general")
import numpy as np
import pandas as pd
from scipy.stats import ConstantInputWarning, spearmanr


def peak_states(values, states):
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all() or len(values) != len(states):
        raise ValueError("Incomplete state support vector")
    return tuple(int(state) for state, value in zip(states, values) if value == values.max())


def support_rho(expression, support):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConstantInputWarning)
        return float(spearmanr(expression, support).statistic)


def summarize_focal_cell_budget(root=REPO, output=None):
    root = Path(root)
    source = root / "results/gse194122"
    output = Path(output) if output is not None else source / "cell_budget"
    supports = pd.read_csv(source / "cell_budget/focal_module_support.tsv", sep="\t")
    expression = pd.read_csv(source / "programs/module_state_summary.tsv.gz", sep="\t")
    heroes = pd.read_csv(source / "programs/hero_regulons.tsv", sep="\t").sort_values("hero_rank")
    draws, summary = [], []
    for module in heroes.module_id:
        rna = expression[expression.module_id == module].sort_values("state_order")
        full = supports[(supports.module_id == module) & (supports.fit_type == "all-cell")].sort_values("state_order")
        if rna.state_order.tolist() != [1, 2, 3, 4, 5] or full.state_order.tolist() != rna.state_order.tolist():
            raise ValueError(f"Incomplete fixed RNA/full-cell reference: {module}")
        stable_peak = peak_states(rna.regulatory_activity, rna.state_order)
        for metric, column in (("target", "selected_targets"), ("chain", "selected_chains")):
            full_values = full[column].to_numpy() / full.module_targets.to_numpy()
            full_peak = peak_states(full_values, full.state_order)
            full_rho = support_rho(rna.median_gene_activity, full_values)
            current = []
            for seed in range(3101, 3111):
                budget = supports[(supports.module_id == module) & (supports.fit_type == "budget-800") & (supports.seed == seed)].sort_values("state_order")
                if budget.state_order.tolist() != rna.state_order.tolist():
                    raise ValueError(f"Incomplete budget support vector: {module}, {seed}")
                values = budget[column].to_numpy() / budget.module_targets.to_numpy()
                peak = peak_states(values, budget.state_order)
                rho = support_rho(rna.median_gene_activity, values)
                row = dict(module_id=module, support_metric=metric, seed=seed,
                    all_cell_peak="|".join(map(str, full_peak)), budget_peak="|".join(map(str, peak)),
                    peak_exactly_preserved=peak == full_peak, peak_overlap=bool(set(peak) & set(full_peak)),
                    budget_peak_tied=len(peak) > 1, all_cell_expression_support_rho=full_rho,
                    budget_expression_support_rho=rho,
                    rho_direction_preserved=bool(np.isfinite(rho) and np.isfinite(full_rho) and np.sign(rho) == np.sign(full_rho)),
                    expression_reference="fixed donor-state RNA module medians from Fig5",
                    support_reference="full-cell fit; no new donor-stable consensus")
                current.append(row)
            draws.extend(current)
            summary.append(dict(module_id=module, support_metric=metric, draws=len(current),
                all_cell_peak="|".join(map(str, full_peak)), donor_stable_chain_peak="|".join(map(str, stable_peak)),
                peak_exact_preserved=sum(x["peak_exactly_preserved"] for x in current),
                peak_overlap=sum(x["peak_overlap"] for x in current),
                tied_budget_peaks=sum(x["budget_peak_tied"] for x in current),
                all_cell_expression_support_rho=full_rho,
                finite_budget_rho=sum(np.isfinite(x["budget_expression_support_rho"]) for x in current),
                rho_direction_preserved=sum(x["rho_direction_preserved"] for x in current)))
    tables = {"focal_budget_draws.tsv": pd.DataFrame(draws), "focal_budget_summary.tsv": pd.DataFrame(summary)}
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".focal-budget-", dir=output) as temporary:
        for name, table in tables.items():
            staged = Path(temporary) / name
            table.to_csv(staged, sep="\t", index=False, lineterminator="\r\n")
            staged.replace(output / name)
    return tables


if __name__ == "__main__":
    summarize_focal_cell_budget()
