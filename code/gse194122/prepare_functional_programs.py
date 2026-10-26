import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "functions"))
from runtime import REPO, use_python

if __name__ == "__main__":
    use_python("general")

import json

import pandas as pd


from erythroid_programs import standardize_donor_expression, reported_program_go


def prepare_functional_programs(root=REPO, output=None):
    root = Path(root).resolve()
    output = Path(output).resolve() if output is not None else root / "results/gse194122/cases"
    output.mkdir(exist_ok=True)
    paths = {
        "expression": "results/gse194122/trajectory/donor_state_expression.tsv.gz",
        "edges": "results/gse194122/donor_stability/stable_consensus_edges.tsv.gz",
        "annotation": "results/gse194122/go_annotation/TERM2GENE.tsv.gz",
        "core": "results/gse194122/cases/HMGA1_gene_recurrence.tsv",
    }
    tables = {key: pd.read_csv(root / path, sep="\t", low_memory=False) for key, path in paths.items()}
    expression = tables["expression"]
    _, standardized = standardize_donor_expression(expression)
    edges, annotation, core = tables["edges"], tables["annotation"], tables["core"]
    definitions = [
        ("KLF1_heme", "KLF1", "GO:0006783", "TF-gene"),
        ("SOX6_positive_erythroid", "SOX6", "GO:0048821", "TF-gene"),
        ("BNC2_chain_immune", "BNC2", "GO:0002252", "factorized-chain"),
        ("HMGA1_recurrent_translation", "HMGA1", "GO:0002181", None),
    ]
    member_rows = []
    for program, tf, term, endpoint in definitions:
        targets = core.loc[core.in_reference & (core.budget_draws >= 8), "target"] if endpoint is None else edges.loc[
            (edges.regulator == tf) & (edges.direction == 1) & (edges.endpoint == endpoint), "target"]
        genes = sorted(set(targets) & set(annotation.loc[annotation.go_id == term, "target"]))
        member_rows.extend(dict(program=program, regulator=tf, direction=1, target=gene) for gene in genes)
    members = pd.DataFrame(member_rows)
    states = {name: i+1 for i, name in enumerate(["HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast"])}
    rows = []
    for program, data in members.groupby("program"):
        tf = data.regulator.iloc[0]
        for signal, scores in [("Target-program RNA", standardized.loc[data.target].mean(axis=0)), ("TF RNA", standardized.loc[tf])]:
            for (state, donor), score in scores.items():
                if pd.notna(score):
                    rows.append(dict(program=program, regulator=tf, genes=len(data), signal=signal, state=state, state_order=states[state], donor_id=donor, RNA_score=score))
    scores = pd.DataFrame(rows)
    scores.to_csv(output / "program_donor_RNA_scores.tsv", sep="\t", index=False)
    members.to_csv(output / "program_members.tsv", sep="\t", index=False)
    genes = sorted(members.loc[members.regulator == "SOX6", "target"])
    def selected_pairs(table):
        selected = table[(table.regulator == "SOX6") & (table.direction == 1) &
                         (table.endpoint == "TF-gene") & table.target.isin(genes)]
        return set(selected[["state_order", "target"]].itertuples(index=False, name=None))
    reference = selected_pairs(edges)
    counts = {(state, gene): 0 for state in range(1, 6) for gene in genes}
    for seed in range(3101, 3111):
        relative = f"results/gse194122/cell_budget/consensus_seed{seed}.tsv.gz"
        for pair in selected_pairs(pd.read_csv(root / relative, sep="\t")):
            counts[pair] += 1
    recurrence = pd.DataFrame([dict(program="SOX6_positive_erythroid", state_order=state, target=gene,
                                   in_reference=(state, gene) in reference, budget_draws=count)
                               for (state, gene), count in sorted(counts.items())])
    recurrence.to_csv(output / "SOX6_signed_link_recurrence.tsv", sep="\t", index=False)
    expression[expression.feature.isin(members.regulator)].to_csv(output / "focal_TF_raw_expression.tsv", sep="\t", index=False)
    metadata = dict(scope="Exploratory GSE194122 program display; no perturbation inputs or network refitting", program_definitions=dict(HMGA1="All 16 cytoplasmic-translation genes among the 23 original loss targets recurring in >=8/10 computational draws; original 23-gene set selected without GO filtering", SOX6="All 14 GO erythrocyte-development genes in the non-self positive direct-target union", KLF1="All 13 GO heme-biosynthesis genes in the non-self positive direct-target union", BNC2="All 39 GO immune-effector-process genes in the non-self positive chain-target union"), RNA_score="Per-feature z scores over the 49 donor-state means; target-program score is their gene mean. TF RNA uses the same per-feature standardization. Both are relative RNA measures from the discovery cohort.")
    reported_program_go(root).to_csv(output / "reported_program_GO.tsv", sep="\t", index=False)
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    prepare_functional_programs()
