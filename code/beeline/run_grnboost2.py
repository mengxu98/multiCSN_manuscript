#!/usr/bin/env python3
"""Run unmodified Arboreto GRNBoost2; no alternative inference backend."""
import json
import os
import sys
from pathlib import Path
import time
from types import SimpleNamespace

# Apply before importing numerical libraries, including in spawned Dask workers.
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code/functions"))
from runtime import use_python
if __name__ == "__main__":
    use_python("grnboost2")
sys.path.insert(0, str(REPO / "code"))
from methods.GRNBoost2 import infer_tf_names, infer_grnboost2



SETTINGS = dict(repo=REPO, output=REPO / "runs/beeline/grnboost2_official", limit=0)


def main():
    import numpy as np
    import pandas as pd
    import arboreto.algo
    import arboreto.core
    from distributed import Client, LocalCluster
    args = SimpleNamespace(**SETTINGS)
    paths = sorted(p for group in ('Curated', 'Synthetic')
                   for p in (args.repo / 'data/beeline' / group).rglob('ExpressionData.csv'))
    assert len(paths) == 500, len(paths)
    if args.limit:
        paths = paths[:args.limit]
    args.output.mkdir(parents=True, exist_ok=False)
    provenance = {
        'backend': 'arboreto.algo.grnboost2', 'seed': 1, 'workers': 1,
        'threads_per_worker': 1,
        'timing_scope': 'Input CSV loading, official GRNBoost2 graph construction/scheduling/inference, and network CSV export; persistent single-worker Dask startup measured separately; imports excluded.',
    }
    (args.output / 'provenance.json').write_text(json.dumps(provenance, indent=2))
    started = time.perf_counter()
    records = []
    with LocalCluster(n_workers=1, threads_per_worker=1, processes=True,
                      dashboard_address=None, memory_limit=0) as cluster, Client(cluster) as client:
        provenance['cluster_startup_seconds'] = time.perf_counter() - started
        (args.output / 'provenance.json').write_text(json.dumps(provenance, indent=2))
        for index, path in enumerate(paths, 1):
            case = path.parent.name
            group = path.parent.parent.name
            typ = path.parent.parent.parent.name.lower()
            dest = args.output / 'networks' / typ / group / (case + '.csv')
            dest.parent.mkdir(parents=True, exist_ok=True)
            start = time.perf_counter()
            expression = pd.read_csv(path, index_col=0).T
            regulators = infer_tf_names(expression, path)
            net = infer_grnboost2(expression, regulators, dest, client, seed=1, verbose=False)
            elapsed = time.perf_counter() - start
            assert len(net) > 0 and np.isfinite(net.weight).all(), case
            assert not net.duplicated(['regulator','target']).any(), case
            assert not (net.regulator == net.target).any(), case
            records.append({'Case':case,'DatasetType':typ,'DatasetGroup':group,
                            'Method':'GRNBoost2','RuntimeSeconds':elapsed,'Status':'complete',
                            'InputCells':len(expression),'InputGenes':len(expression.columns),
                            'InputRegulators':len(regulators),'OutputEdges':len(net),
                            'TimingScope':provenance['timing_scope']})
            pd.DataFrame(records).to_csv(args.output / 'runtime.tsv', sep='\t', index=False)
            (args.output / 'progress.json').write_text(json.dumps({'completed':index,'expected':len(paths),'last_case':case,'total_case_seconds':sum(r['RuntimeSeconds'] for r in records),'load':os.getloadavg()}))
            print(f'{index}/{len(paths)} {case} {elapsed:.3f}s {len(net)} edges', flush=True)
    (args.output / 'complete.json').write_text(json.dumps({'cases':len(records),'wall_seconds':time.perf_counter()-started}))


if __name__ == '__main__':
    main()
