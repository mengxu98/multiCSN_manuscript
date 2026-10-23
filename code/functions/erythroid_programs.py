from pathlib import Path

import pandas as pd


def standardize_donor_expression(expression, columns=None):
    wide = expression.pivot(index="feature", columns=["state", "donor_id"], values="mean_normalized_expression")
    if columns is not None:
        wide = wide.reindex(columns=columns)
    standardized = wide.sub(wide.mean(axis=1), axis=0).div(wide.std(axis=1, ddof=1).replace(0, 1), axis=0)
    return wide, standardized


def reported_program_go(root):
    root = Path(root)
    records = []
    for endpoint, tag, module, term in [
        ("TF-gene", "direct", "SOX6 (+)", "GO:0048821"),
        ("TF-gene", "direct", "KLF1 (+)", "GO:0006783"),
        ("factorized-chain", "chain", "BNC2 (+)", "GO:0002252"),
    ]:
        table = pd.read_csv(root / f"results/gse194122/programs/{tag}_module_GO_supported.tsv", sep="\t")
        selected = table[(table.module_id == module) & (table.go_id == term)]
        record = selected.iloc[0].to_dict()
        record["p_adjust"] = record.pop("q_global")
        record["endpoint"] = endpoint
        record["test_family"] = "Global BH across all signed programs and eligible GO terms within endpoint"
        records.append(record)
    table = pd.read_csv(root / "results/gse194122/cases/HMGA1_GO_all_tests.tsv.gz", sep="\t")
    selected = table[(table["set"] == "core_8_of_10") & (table.background == "HMGA1_source_targets") & (table.go_id == "GO:0002181")]
    record = selected.iloc[0].to_dict()
    record.update(module_id="HMGA1 recurrent loss core (23 genes)",
                  annotated_set=record["annotated_set_genes"],
                  annotated_universe=record["annotated_background_genes"],
                  endpoint="TF-gene loss event",
                  test_family="BH across all eligible terms for the 23-gene core against fixed HMGA1 source-state targets; RNA plot uses its 16 translation members")
    records.append(record)
    return pd.DataFrame(records)[["module_id", "go_id", "term", "overlap", "term_genes", "annotated_set", "annotated_universe", "p_value", "p_adjust", "endpoint", "test_family"]]
