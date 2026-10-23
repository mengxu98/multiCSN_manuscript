#!/usr/bin/env python3
"""Run the BEELINE benchmark with one numerical thread per method."""
import os
from pathlib import Path
import random
import subprocess
import sys
import time
from types import SimpleNamespace

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
os.environ["LC_ALL"] = "C"

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code/functions"))
from runtime import rscript_executable, PYTHON_LOCATIONS, use_python, r_function_command, python_function_command
METHODS = ("multiCSN", "PPCOR", "SCORPION", "PIDC", "LEAP", "GRNBoost2",
           "GNIPLR", "SCODE", "SCSGL", "SINCERITIES", "GENIE3")
PY_METHODS = {"GNIPLR", "GRNBoost2", "SCSGL"}


def case_paths(limit):
    # Match R's bytewise full-path order (Path's component ordering differs for
    # prefix case names such as case-1 versus case-1-10).
    paths = sorted((p for group in ("Curated", "Synthetic")
                    for p in (REPO / "data/beeline" / group).rglob("ExpressionData.csv")), key=str)
    if len(paths) != 500:
        raise ValueError(f"Expected 500 inputs, got {len(paths)}")
    return paths[:limit] if limit else paths


def python_worker(output, method, limit):
    output = Path(output)
    import numpy as np
    import pandas as pd
    from contextlib import ExitStack
    sys.path.insert(0, str(REPO / "code"))
    from methods.GNIPLR import infer_gniplr
    from methods.SCSGL import infer_scsgl
    runners = {"GNIPLR": infer_gniplr, "SCSGL": infer_scsgl}
    records = []
    with ExitStack() as stack:
        if method == "GRNBoost2":
            from arboreto.algo import grnboost2
            from distributed import Client, LocalCluster
            from methods.GRNBoost2 import infer_tf_names, infer_grnboost2
            cluster = stack.enter_context(LocalCluster(n_workers=1, threads_per_worker=1,
                                                       processes=True, dashboard_address=None,
                                                       memory_limit=0))
            client = stack.enter_context(Client(cluster))
        # Import/warm up code, but do not run a case outside the measured boundary.
        if method == "GNIPLR":
            import scipy.stats
            import sklearn.linear_model
            import statsmodels.api
        paths = case_paths(limit)
        for index, path in enumerate(paths, 1):
            case, group, typ = path.parent.name, path.parent.parent.name, path.parent.parent.parent.name.lower()
            relative = Path("networks") / f"{method}__{case}.csv"
            dest = output / relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            random.seed(1)
            np.random.seed(1)
            # Reading the case and preparing the in-memory matrix are excluded from
            # the measured boundary; only inference and network export are timed.
            load_started = time.perf_counter()
            if method == "GRNBoost2":
                expression = pd.read_csv(path, index_col=0).T
                regulators = infer_tf_names(expression, path)
            else:
                expression = pd.read_csv(path, index_col=0)
            load_seconds = time.perf_counter() - load_started
            started = time.perf_counter()
            if method == "GRNBoost2":
                infer_grnboost2(expression, regulators, dest, client, seed=1, verbose=False)
            else:
                runners[method](expression, dest)
            seconds = time.perf_counter() - started
            records.append(dict(Case=case, Method=method, DatasetType=typ, DatasetGroup=group,
                                RuntimeSeconds=seconds, LoadSeconds=load_seconds,
                                LoadAverage=os.getloadavg()[0],
                                Network=str(relative),
                                TimingScope="inference_and_export", Seed=1, Status="complete"))
            pd.DataFrame(records).to_csv(output / "work" / f"{method}.tsv", sep="\t", index=False)
            print(f"{index}/{len(paths)} {method} {case} {seconds:.3f}s", flush=True)


def validate(output, methods, limit):
    import numpy as np
    import pandas as pd
    frames = []
    expected = {(p.parent.name, p.parent.parent.name, p.parent.parent.parent.name.lower())
                for p in case_paths(limit)}
    for method in methods:
        table = pd.read_csv(output / "work" / f"{method}.tsv", sep="\t")
        actual = set(zip(table.Case, table.DatasetGroup, table.DatasetType))
        if len(table) != len(expected) or actual != expected or set(table.Method) != {method}:
            raise ValueError(f"Incomplete case contract: {method}")
        if not np.isfinite(table.RuntimeSeconds).all() or (table.RuntimeSeconds <= 0).any():
            raise ValueError(f"Invalid timings: {method}")
        if not table.Status.eq("complete").all():
            raise ValueError(f"Unfinished records: {method}")
        missing_weights = []
        for row in table.itertuples():
            expected_path = Path("networks") / f"{method}__{row.Case}.csv"
            if Path(row.Network) != expected_path:
                raise ValueError(f"Case/network identity mismatch: {row.Case}")
            path = output / row.Network
            net = pd.read_csv(path).rename(columns={"regulatoryGene": "regulator", "targetGene": "target"})
            if list(net.columns) != ["regulator", "target", "weight"] or net.empty:
                raise ValueError(f"Invalid network: {path}")
            if net[["regulator", "target"]].isna().any().any() or net.duplicated(["regulator", "target"]).any():
                raise ValueError(f"Invalid endpoints: {path}")
            # Preserve method-defined undefined weights, reporting rather than silently replacing them.
            weights = pd.to_numeric(net.weight, errors="raise")
            if np.isinf(weights).any():
                raise ValueError(f"Infinite weights: {path}")
            missing_weights.append(int(weights.isna().sum()))
        table["UndefinedWeights"] = missing_weights
        frames.append(table)
    combined = pd.concat(frames, ignore_index=True)
    (output / "runtime").mkdir(exist_ok=True)
    runtime_columns = ["Case", "Method", "DatasetType", "DatasetGroup", "RuntimeSeconds"]
    if "LoadSeconds" in combined.columns:
        runtime_columns.append("LoadSeconds")
    if "LoadAverage" in combined.columns:
        runtime_columns.append("LoadAverage")
    runtime_columns.append("TimingScope")
    combined[runtime_columns].to_csv(output / "runtime/case_runtime.tsv", sep="\t", index=False)
    summary = combined.groupby("Method").RuntimeSeconds.agg(
        cases="size", total_runtime_seconds="sum", median_runtime_seconds="median", maximum_runtime_seconds="max")
    summary["runtime_ratio_to_multicsn"] = summary.total_runtime_seconds / summary.loc["multiCSN", "total_runtime_seconds"]
    summary.sort_values("total_runtime_seconds").to_csv(output / "runtime/runtime_summary.tsv", sep="\t")
    return len(combined)


SETTINGS = dict(output=REPO / "runs/beeline/unified_formal", limit=0)


def main():
    use_python("beeline")
    args = SimpleNamespace(**SETTINGS, rscript=rscript_executable("beeline"),
                           grnboost2_python=str(PYTHON_LOCATIONS["grnboost2"][0]))
    if args.limit < 0:
        raise ValueError("limit must be nonnegative")
    output = args.output.resolve()
    os.chdir(REPO)
    if output == REPO / "results/beeline" or REPO / "results" in output.parents:
        raise ValueError("Use an isolated staging directory outside formal results")
    output.mkdir(parents=True, exist_ok=False)
    (output / "work").mkdir()
    os.environ["PATH"] = str(Path(args.rscript).resolve().parent) + os.pathsep + os.environ["PATH"]
    for method in METHODS:

        if method in PY_METHODS:
            executable = args.grnboost2_python if method == "GRNBoost2" else sys.executable
            command = python_function_command(
                "code/beeline/run_unified_benchmark.py", "python_worker", executable=executable,
                output=str(output), method=method, limit=args.limit)
        else:
            command = r_function_command(
                "code/beeline/unified_worker.R", "run_benchmark_method", executable=args.rscript,
                output_root=output, method=method, limit=args.limit)
        subprocess.run(command, check=True)
    validate(output, METHODS, args.limit)
    if not args.limit:

        for script in ("reevaluate_networks.R", "prepare_complete_method_tables.R"):
            subprocess.run(r_function_command(
                "code/beeline/" + script, script.removesuffix(".R"), executable=args.rscript,
                result_root=output), check=True)



if __name__ == "__main__":
    main()
