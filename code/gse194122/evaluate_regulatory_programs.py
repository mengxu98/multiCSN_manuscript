import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "functions"))
from runtime import REPO, use_python

if __name__ == "__main__":
    use_python("general")


import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import hypergeom, wilcoxon
from erythroid_programs import standardize_donor_expression


def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind="stable")
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(adjusted, 1)
    return result


def evaluate_regulatory_programs(root=REPO, output=None):
    root = Path(root).resolve()
    out = Path(output).resolve() if output is not None else root / "results/gse194122/programs"
    out.mkdir(parents=True, exist_ok=True)


    def read(relative):
        path = root / relative
        return pd.read_csv(path, sep="\t", low_memory=False)


    edges = read("results/gse194122/donor_stability/stable_consensus_edges.tsv.gz")
    expression = read("results/gse194122/trajectory/donor_state_expression.tsv.gz")
    groups = read("results/gse194122/trajectory/donor_state_group_counts.tsv")
    states = ["HSC", "MK/E prog", "Proerythroblast", "Erythroblast", "Normoblast"]
    go = read("results/gse194122/go_annotation/TERM2GENE.tsv.gz").drop_duplicates()
    names = read("results/gse194122/go_annotation/TERM2NAME.tsv").set_index("go_id").term
    primary = read("results/gse194122/go_annotation/displayed_primary_terms.tsv")
    groups = groups.assign(state_order=groups.state.map(dict(zip(states, range(1, 6))))).sort_values(["state_order", "donor_id"])
    columns = pd.MultiIndex.from_frame(groups[["state", "donor_id"]])
    wide, z = standardize_donor_expression(expression, columns)
    universe = set(edges.loc[edges.endpoint == "TF-gene", "target"])
    mapping = go[go.target.isin(universe)]
    sizes = mapping.groupby("go_id").target.nunique()
    terms = sizes[sizes.between(10, 500)].index
    mapping = mapping[mapping.go_id.isin(terms)]
    genes = sorted(set(mapping.target))
    gi = {g: i for i, g in enumerate(genes)}
    ti = {t: i for i, t in enumerate(terms)}
    term_matrix = csr_matrix((np.ones(len(mapping), dtype=np.int32), ([ti[t] for t in mapping.go_id], [gi[g] for g in mapping.target])), shape=(len(terms), len(genes)))
    term_sizes = np.asarray(term_matrix.sum(axis=1)).ravel()
    module_rows = []
    score_rows = []
    contrast_rows = []
    support_rows = []
    for endpoint in ["TF-gene", "factorized-chain"]:
        e = edges[(edges.endpoint == endpoint) & (edges.regulator != edges.target)].copy()
        e["module_id"] = e.regulator + np.where(e.direction > 0, " (+)", " (-)")
        memberships = e.groupby("module_id").target.agg(lambda x: sorted(set(x)))
        memberships = memberships[memberships.map(len) >= 10]
        annotated = csr_matrix((np.ones(sum(len(set(v) & set(gi)) for v in memberships), dtype=np.int32),
                               ([i for i, v in enumerate(memberships) for g in v if g in gi], [gi[g] for v in memberships for g in v if g in gi])),
                              shape=(len(memberships), len(genes)))
        overlaps = (annotated @ term_matrix.T).toarray()
        n = np.asarray(annotated.sum(axis=1)).ravel()
        p = hypergeom.sf(overlaps - 1, len(genes), term_sizes[None, :], n[:, None])
        q = bh(p.ravel()).reshape(p.shape)
        go_parts = []
        for i, (module, targets) in enumerate(memberships.items()):
            score = z.loc[targets].mean(axis=0)
            state_score = score.groupby(level="state").median().reindex(states)
            donor_score = score.unstack("state").reindex(columns=states)
            peak = state_score.idxmax()
            peak_count = int((donor_score.idxmax(axis=1) == peak).sum())
            tf = module.rsplit(" (", 1)[0]
            module_rows.append(dict(endpoint=endpoint, module_id=module, regulator=tf, targets=len(targets), RNA_peak=peak,
                                    peak_donors=peak_count, available_donors=len(donor_score), RNA_peak_value=state_score.max(), RNA_trough_value=state_score.min()))
            for (state, donor), value in score.items():
                score_rows.append(dict(endpoint=endpoint, module_id=module, state=state, donor_id=donor, score=value))
            for j in range(4):
                delta = (donor_score[states[j+1]] - donor_score[states[j]]).dropna()
                pvalue = wilcoxon(delta, alternative="two-sided", method="auto").pvalue if (delta != 0).any() else 1
                contrast_rows.append(dict(endpoint=endpoint, module_id=module, transition_order=j+1, donors=len(delta),
                                          median_change=delta.median(), positive_donors=int((delta > 0).sum()), negative_donors=int((delta < 0).sum()), p_value=pvalue))
            gp = pd.DataFrame(dict(module_id=module, go_id=terms, term=[names.get(t, t) for t in terms], overlap=overlaps[i],
                                   term_genes=term_sizes, annotated_set=n[i], annotated_universe=len(genes), p_value=p[i], q_global=q[i]))
            go_parts.append(gp)
            for state in states:
                paths = e[(e.module_id == module) & (e.state == state)]
                support_rows.append(dict(endpoint=endpoint, module_id=module, state=state, paths=len(paths), targets=paths.target.nunique(), regions=paths.region.nunique()))
        go_tests = pd.concat(go_parts, ignore_index=True)
        tag = "direct" if endpoint == "TF-gene" else "chain"
        go_tests.to_csv(out / f"{tag}_module_GO_all_tests.tsv.gz", sep="\t", index=False)
        go_tests[go_tests.q_global < 0.05].to_csv(out / f"{tag}_module_GO_supported.tsv", sep="\t", index=False)
        print(endpoint, "modules", len(memberships), "significant module-term pairs", int((q < .05).sum()), flush=True)
    pd.DataFrame(module_rows).to_csv(out / "all_module_RNA_summary.tsv", sep="\t", index=False)
    pd.DataFrame(score_rows).to_csv(out / "all_module_donor_scores.tsv.gz", sep="\t", index=False)
    contrasts = pd.DataFrame(contrast_rows)
    contrasts["q_global"] = contrasts.groupby("endpoint").p_value.transform(bh)
    contrasts.to_csv(out / "all_module_RNA_contrasts.tsv", sep="\t", index=False)
    pd.DataFrame(support_rows).to_csv(out / "all_module_support.tsv", sep="\t", index=False)
    case_rows = []
    case_terms = []
    for tf, transition, event in [("HMGA1", 4, "lost"), ("BNC2", 3, "gained"), ("BNC2", 4, "lost")]:
        e = edges[(edges.endpoint == "TF-gene") & (edges.regulator == tf)]
        a = set(e.loc[e.state_order == transition, "target"])
        b = set(e.loc[e.state_order == transition + 1, "target"])
        selected = a - b if event == "lost" else b - a
        own = a if event == "lost" else b
        case_rows.append(dict(regulator=tf, transition_order=transition, event=event, previous_targets=len(a), next_targets=len(b), event_targets=len(selected), fraction_own_state=len(selected)/len(own), targets=";".join(sorted(selected))))
        for background, bg in [("full_stable_targets", universe), ("own_source_or_destination", own)]:
            m = go[go.target.isin(bg)]
            sz = m.groupby("go_id").target.nunique()
            sz = sz[sz.between(10, 500)]
            m = m[m.go_id.isin(sz.index)]
            eligible = set(m.target)
            hits = m[m.target.isin(selected)].groupby("go_id").target.nunique().reindex(sz.index, fill_value=0)
            pv = hypergeom.sf(hits.to_numpy()-1, len(eligible), sz.to_numpy(), len(selected & eligible))
            frame = pd.DataFrame(dict(regulator=tf, transition_order=transition, event=event, background=background, go_id=sz.index,
                                      term=[names.get(t,t) for t in sz.index], overlap=hits.to_numpy(), term_genes=sz.to_numpy(), annotated_set=len(selected & eligible), annotated_background=len(eligible), p_value=pv, q_family=bh(pv)))
            case_terms.append(frame)
    pd.DataFrame(case_rows).to_csv(out / "HMGA1_BNC2_event_sizes.tsv", sep="\t", index=False)
    case_go = pd.concat(case_terms, ignore_index=True)
    case_go.to_csv(out / "HMGA1_BNC2_background_GO_all_tests.tsv.gz", sep="\t", index=False)
    case_go.merge(primary[["regulator", "transition_order", "go_id"]].drop_duplicates(), on=["regulator", "transition_order", "go_id"]).to_csv(out / "HMGA1_BNC2_original_terms.tsv", sep="\t", index=False)
    tf_expression = expression[expression.feature.isin(edges.regulator.unique())].copy()
    tf_expression.to_csv(out / "all_TF_donor_expression.tsv.gz", sep="\t", index=False)
    print("PASS: exploratory program tables generated; inputs unchanged", flush=True)


if __name__ == "__main__":
    evaluate_regulatory_programs()
