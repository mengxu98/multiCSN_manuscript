import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "functions"))
from runtime import REPO, use_python

if __name__ == "__main__":
    use_python("general")


import pandas as pd



def summarize_link_recurrence(root=REPO, output=None):
    root = Path(root).resolve()
    out = Path(output).resolve() if output is not None else root / "results/gse194122/cases"
    out.mkdir(exist_ok=True)


    def read(path):
        return pd.read_csv(path, sep="\t", low_memory=False)


    ref = read(root / "results/gse194122/donor_stability/stable_consensus_edges.tsv.gz")
    draws = {"reference": ref}
    for seed in range(3101, 3111):
        draws[str(seed)] = read(root / f"results/gse194122/cell_budget/consensus_seed{seed}.tsv.gz")
    rows = []
    events = []
    counts = []
    focal = ["HMGA1", "BNC2", "KLF1", "NFE2", "SOX6", "GATA1", "GATA2", "MYB", "MECOM", "GFI1B", "SPI1"]
    for draw, d in draws.items():
        for endpoint in ["TF-gene", "factorized-chain"]:
            x = d[d.endpoint == endpoint]
            totals = x.groupby("state_order").size()
            c = x.groupby(["regulator", "state_order"]).agg(links=("target", "size"), targets=("target", "nunique")).reset_index()
            regulators = sorted(set(ref.loc[ref.endpoint == endpoint, "regulator"]))
            grid = pd.MultiIndex.from_product([regulators, range(1, 6)], names=["regulator", "state_order"])
            c = c.set_index(["regulator", "state_order"]).reindex(grid, fill_value=0).reset_index()
            c["network_fraction"] = c.links / c.state_order.map(totals)
            c["draw"] = draw
            c["endpoint"] = endpoint
            counts.append(c)
            for tf in focal:
                own = x[x.regulator == tf]
                for t in range(1, 5):
                    previous = set(own.loc[own.state_order == t, "target"])
                    following = set(own.loc[own.state_order == t+1, "target"])
                    rows.append(dict(draw=draw, endpoint=endpoint, regulator=tf, transition_order=t,
                                     previous_targets=len(previous), next_targets=len(following), gained=len(following-previous),
                                     lost=len(previous-following), retained=len(previous&following)))
                    if endpoint == "TF-gene":
                        for event, members in [("gained", following-previous), ("lost", previous-following)]:
                            for gene in sorted(members):
                                events.append(dict(draw=draw, regulator=tf, transition_order=t, event=event, target=gene))
    pd.concat(counts, ignore_index=True).to_csv(out / "all_TF_state_support.tsv.gz", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(out / "focal_TF_transition_sizes.tsv", sep="\t", index=False)
    ev = pd.DataFrame(events)
    keys = ["regulator", "transition_order", "event", "target"]
    rec = ev[ev.draw != "reference"].groupby(keys).draw.nunique().rename("budget_draws").reset_index()
    rec = ev[ev.draw == "reference"][keys].merge(rec, on=keys, how="left").fillna({"budget_draws":0})
    rec.to_csv(out / "focal_reference_gene_event_recurrence.tsv.gz", sep="\t", index=False)
    expression = read(root / "results/gse194122/trajectory/donor_state_expression.tsv.gz")
    states = ["HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast"]
    own_expression = expression[expression.feature.isin(ref.regulator.unique())]
    stats = own_expression.groupby(["feature", "state"]).agg(donors=("donor_id", "nunique"), median_expression=("mean_normalized_expression", "median"),
        median_detection=("detection_fraction", "median"), nonzero_donors=("mean_normalized_expression", lambda x:int((x>0).sum()))).reset_index()
    stats.to_csv(out / "all_TF_expression_by_state.tsv", sep="\t", index=False)
    contrasts = []
    for tf, data in own_expression.groupby("feature"):
        matrix = data.pivot(index="donor_id", columns="state", values="mean_normalized_expression").reindex(columns=states)
        for t in range(4):
            delta = (matrix[states[t+1]]-matrix[states[t]]).dropna()
            contrasts.append(dict(regulator=tf, transition_order=t+1, donors=len(delta), median_change=delta.median(), positive_donors=int((delta>0).sum()),
                                  negative_donors=int((delta<0).sum()), zero_donors=int((delta==0).sum())))
    pd.DataFrame(contrasts).to_csv(out / "all_TF_paired_expression_changes.tsv", sep="\t", index=False)
    chain = ref[ref.endpoint == "factorized-chain"]
    direct = ref[ref.endpoint == "TF-gene"]
    switches = []
    for t in range(1, 5):
        a = chain[chain.state_order == t]
        b = chain[chain.state_order == t+1]
        pairs = a.merge(b, on=["regulator", "target", "direction"], suffixes=("_previous", "_next"))
        pairs = pairs[pairs.region_previous != pairs.region_next].copy()
        for row in pairs.itertuples():
            before = direct[(direct.state_order == t)&(direct.regulator==row.regulator)&(direct.target==row.target)]
            after = direct[(direct.state_order == t+1)&(direct.regulator==row.regulator)&(direct.target==row.target)]
            switches.append(dict(regulator=row.regulator, target=row.target, direction=row.direction, transition_order=t,
                                 region_previous=row.region_previous, region_next=row.region_next,
                                 direct_pair_present_both=len(before)>0 and len(after)>0,
                                 direct_pair_same_sign=bool(len(before)>0 and len(after)>0 and before.direction.iloc[0]==after.direction.iloc[0])))
    pd.DataFrame(switches).to_csv(out / "full_cohort_candidate_regional_changes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    summarize_link_recurrence()
