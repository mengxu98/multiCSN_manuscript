#!/usr/bin/env python3
"""Stage one auditable BEELINE batch. Never replaces formal results automatically.

The task selects its declared Python and R environments.
All methods run sequentially with one numerical thread; The declared limit is for smoke tests only.
"""
import csv
import importlib.metadata
import json
import os
from pathlib import Path
import platform
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
from runtime import rscript_executable, PYTHON_LOCATIONS, use_python, sha256, r_function_command, python_function_command
METHODS = ("multiCSN", "PPCOR", "SCORPION", "PIDC", "LEAP", "GRNBoost2",
           "GNIPLR", "SCODE", "SCSGL", "SINCERITIES", "GENIE3")
PY_METHODS = {"GNIPLR", "GRNBoost2", "SCSGL"}


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True))
    temporary.replace(path)


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
    write_json(output / "work" / f"{method}_environment.json", dict(
        executable=sys.executable, python=sys.version,
        packages={d.metadata["Name"]: d.version for d in importlib.metadata.distributions()}))
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
                                Network=str(relative), NetworkSHA256=sha256(dest),
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
            if sha256(path) != row.NetworkSHA256:
                raise ValueError(f"Network hash mismatch: {path}")
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
    combined.to_csv(output / "provenance/run_manifest.tsv", sep="\t", index=False)
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
    (output / "provenance").mkdir()
    os.environ["PATH"] = str(Path(args.rscript).resolve().parent) + os.pathsep + os.environ["PATH"]
    inputs = sorted({p for expr in case_paths(args.limit) for p in expr.parent.glob("*.csv")})
    sources = sorted(p for p in (REPO / "code").rglob("*")
                     if p.is_file() and p.suffix in {".R", ".r", ".py", ".sh", ".cpp", ".h"})
    fingerprints = {str(p.relative_to(REPO)): sha256(p) for p in inputs + sources}
    write_json(output / "provenance/input_code_sha256.json", fingerprints)
    write_json(output / "provenance/environment.json", dict(
        python=sys.version, executable=sys.executable, rscript=args.rscript,
        grnboost2_python=args.grnboost2_python,
        platform=platform.platform(), cpu=platform.processor(), cpu_count=os.cpu_count(),
        affinity=sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
        packages={d.metadata["Name"]: d.version for d in importlib.metadata.distributions()},
        load_at_start=os.getloadavg(), seed=1, threads=1, sequential_methods=True,
        timing_scope=("Input reading and matrix preparation are excluded from the measured region; "
                      "timing covers inference and network export only (TimingScope=inference_and_export). "
                      "Per-case LoadSeconds and LoadAverage are recorded. Process startup and persistent "
                      "Dask startup excluded; lazy imports inside method wrappers included."),
        seed_limitations="SCODE external R subprocess currently uses its own RNG; parent seed does not control it.",
        settings={"output": str(output), "limit": args.limit}, started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())))
    complete = []
    try:
        for method in METHODS:
            write_json(output / "status.json", dict(state="running", method=method, completed_methods=complete))
            if method in PY_METHODS:
                executable = args.grnboost2_python if method == "GRNBoost2" else sys.executable
                command = python_function_command(
                    "code/beeline/run_unified_benchmark.py", "python_worker", executable=executable,
                    output=str(output), method=method, limit=args.limit)
            else:
                command = r_function_command(
                    "code/beeline/unified_worker.R", "run_benchmark_method", executable=args.rscript,
                    output_root=output, method=method, limit=args.limit)
            with (output / "work" / f"{method}.log").open("w") as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
            complete.append(method)
        changed = [name for name, digest in fingerprints.items() if sha256(REPO / name) != digest]
        if changed:
            raise ValueError(f"Inputs/code changed during run: {changed[:10]}")
        count = validate(output, METHODS, args.limit)
        environment_path = output / "provenance/environment.json"
        environment = json.loads(environment_path.read_text())
        environment["python_workers"] = {
            method: json.loads((output / "work" / f"{method}_environment.json").read_text())
            for method in PY_METHODS}
        environment["r_sessions"] = {
            method: (output / "work" / f"{method}_session.txt").read_text()
            for method in METHODS if method not in PY_METHODS}
        with (output / "work/multiCSN_packages.tsv").open() as stream:
            environment["r_packages"] = list(csv.DictReader(stream, delimiter="\t"))
        environment["load_at_end"] = os.getloadavg()
        write_json(environment_path, environment)
        if not args.limit:
            write_json(output / "status.json", dict(state="evaluating", cases=count, completed_methods=complete))
            for script in ("reevaluate_networks.R", "prepare_complete_method_tables.R"):
                with (output / "work" / (script + ".log")).open("w") as log:
                    subprocess.run(r_function_command(
                        "code/beeline/" + script, script.removesuffix(".R"), executable=args.rscript,
                        result_root=output), stdout=log, stderr=subprocess.STDOUT, check=True)
        write_json(output / "status.json", dict(state="validated", cases=count,
                   smoke_test=bool(args.limit), completed_methods=complete, promoted=False))
    except Exception as exc:
        write_json(output / "status.json", dict(state="failed", error=str(exc), completed_methods=complete))
        raise


if __name__ == "__main__":
    main()
