#!/usr/bin/env python3
"""Validate an isolated unified BEELINE batch and install it into results/beeline.

Defaults are declared in this script. Existing formal results are backed up
before admission; core_scaling remains an independent experiment.
"""
import json
import math
import shutil
import tempfile
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "functions"))
from runtime import use_python, sha256
from file_io import publish_directory
sys.path.insert(0, str(Path(__file__).absolute().parent))
from run_unified_benchmark import METHODS, case_paths
if __name__ == "__main__":
    use_python("beeline")

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
STAGE = REPO / "runs/beeline/unified_formal"
BACKUP_ROOT = Path.home() / "backups/multiCSN_manuscript/beeline_promotion"
TARGET = REPO / "results" / "beeline"
KEEP_IN_PLACE = {"core_scaling"}
EXPECTED_METHODS = len(METHODS)
EXPECTED_CASES = 500
EXPECTED_RECORDS = EXPECTED_METHODS * EXPECTED_CASES
METRIC_FILES = (
    "curated__primary_metrics.csv", "curated__other_metrics.csv",
    "synthetic__primary_metrics.csv", "synthetic__other_metrics.csv",
    "combined__primary_metrics.csv", "combined__other_metrics.csv",
    "case_metrics.tsv", "case_ranks.tsv", "rank_summary.tsv",
)


def fail(message):
    raise SystemExit(f"promotion failed: {message}")


def validate_batch():
    required = ["provenance/run_manifest.tsv", "provenance/environment.json",
                "provenance/input_code_sha256.json", "runtime/case_runtime.tsv",
                "runtime/runtime_summary.tsv", "status.json",
                "evaluation/case_metrics.tsv", "evaluation/source_index.tsv",
                "evaluation/metric_policy.tsv", "evaluation/beeline_method_summary.tsv",
                "evaluation/beeline_case_manifest.tsv"]
    required += [f"metrics/{name}" for name in METRIC_FILES]
    for name in required:
        if not (STAGE / name).is_file():
            fail(f"missing required batch file: {name}")
    hashes = {name: sha256(STAGE / name) for name in required}
    status_path = STAGE / "status.json"
    status = json.loads(status_path.read_text())
    if status.get("state") != "validated":
        fail(f"batch state is not 'validated': {status}")
    if status.get("smoke_test"):
        fail("smoke-test batches cannot be promoted")
    if status.get("cases") != EXPECTED_RECORDS or set(status.get("completed_methods", [])) != set(METHODS):
        fail("batch status does not cover the declared methods and cases")

    manifest = pd.read_csv(STAGE / "provenance" / "run_manifest.tsv", sep="\t")
    if len(manifest) != EXPECTED_RECORDS:
        fail(f"run manifest has {len(manifest)} rows, expected {EXPECTED_RECORDS}")
    runtime = pd.read_csv(STAGE / "runtime" / "case_runtime.tsv", sep="\t")
    if len(runtime) != EXPECTED_RECORDS:
        fail(f"case runtime has {len(runtime)} rows, expected {EXPECTED_RECORDS}")
    columns = ["Case", "Method", "DatasetType", "DatasetGroup", "RuntimeSeconds", "TimingScope"]
    registered = {(path.parent.name, path.parent.parent.parent.name.lower(), path.parent.parent.name)
                  for path in case_paths(0)}
    if len(registered) != EXPECTED_CASES:
        fail("registered case inventory is incomplete")
    for label, table in (("manifest", manifest), ("runtime", runtime)):
        if not set(columns).issubset(table.columns) or table[columns].isna().any().any():
            fail(f"invalid {label} schema or missing values")
        if table.duplicated(["Method", "Case"]).any() or set(table.Method) != set(METHODS):
            fail(f"duplicate or unexpected {label} case/method identities")
        cases = [set(zip(group.Case, group.DatasetType, group.DatasetGroup))
                 for _, group in table.groupby("Method")]
        if any(group != registered for group in cases):
            fail(f"incomplete {label} case coverage")
        if not table.RuntimeSeconds.map(lambda value: math.isfinite(value) and value > 0).all():
            fail(f"invalid {label} runtimes")
        if not table.TimingScope.eq("inference_and_export").all():
            fail(f"invalid {label} timing scope")
    if not {"Status", "Network", "NetworkSHA256"}.issubset(manifest.columns) or not manifest.Status.eq("complete").all():
        fail("incomplete network manifest")
    for optional in ("LoadSeconds", "LoadAverage"):
        if optional in manifest.columns or optional in runtime.columns:
            columns.append(optional)
    try:
        pd.testing.assert_frame_equal(
            manifest.sort_values(["Method", "Case"])[columns].reset_index(drop=True),
            runtime.sort_values(["Method", "Case"])[columns].reset_index(drop=True),
            check_dtype=False, check_exact=True)
    except (AssertionError, KeyError) as error:
        fail(f"runtime and network manifest disagree: {error}")

    summary = pd.read_csv(STAGE / "runtime/runtime_summary.tsv", sep="\t").set_index("Method")
    expected = runtime.groupby("Method").RuntimeSeconds.agg(
        cases="size", total_runtime_seconds="sum", median_runtime_seconds="median", maximum_runtime_seconds="max")
    expected["runtime_ratio_to_multicsn"] = expected.total_runtime_seconds / expected.loc["multiCSN", "total_runtime_seconds"]
    try:
        pd.testing.assert_frame_equal(summary.sort_index(), expected.sort_index(),
                                      check_dtype=False, rtol=1e-12, atol=1e-12)
    except AssertionError as error:
        fail(f"runtime summary disagrees with case runtimes: {error}")
    for name in ("environment.json", "input_code_sha256.json"):
        value = json.loads((STAGE / "provenance" / name).read_text())
        if not isinstance(value, dict) or not value:
            fail(f"empty or invalid provenance: {name}")
    evaluated = pd.read_csv(STAGE / "evaluation/case_metrics.tsv", sep="\t")
    primary = ["AUPRC", "EarlyPrecision", "EPR"]
    identity = ["Case", "Method", "DatasetType", "DatasetGroup"]
    try:
        pd.testing.assert_frame_equal(
            evaluated.sort_values(identity)[identity].reset_index(drop=True),
            manifest.sort_values(identity)[identity].reset_index(drop=True),
            check_dtype=False, check_exact=True)
        if not evaluated[primary].apply(lambda column: column.map(math.isfinite)).all().all():
            fail("nonfinite evaluated primary metrics")
        for name in ("case_metrics.tsv", "combined__primary_metrics.csv"):
            table = pd.read_csv(STAGE / "metrics" / name,
                                sep="\t" if name.endswith(".tsv") else ",")
            columns = identity + primary
            pd.testing.assert_frame_equal(
                table.sort_values(identity)[columns].reset_index(drop=True),
                evaluated.sort_values(identity)[columns].reset_index(drop=True),
                check_dtype=False, rtol=1e-12, atol=1e-12)
        ranks = pd.read_csv(STAGE / "metrics/case_ranks.tsv", sep="\t")
        expected_ranks = evaluated.melt(id_vars=identity, value_vars=primary,
                                       var_name="Metric", value_name="Value")
        expected_ranks["Rank"] = expected_ranks.groupby(["Case", "Metric"]).Value.rank(
            ascending=False, method="average")
        keys = identity + ["Metric"]
        columns = keys + ["Value", "Rank"]
        pd.testing.assert_frame_equal(
            ranks.sort_values(keys)[columns].reset_index(drop=True),
            expected_ranks.sort_values(keys)[columns].reset_index(drop=True),
            check_dtype=False, rtol=1e-12, atol=1e-12)
    except (AssertionError, KeyError, TypeError) as error:
        fail(f"evaluation and plotting tables disagree: {error}")
    for row in manifest.itertuples():
        name = f"networks/{row.Method}__{row.Case}.csv"
        if row.Network != name or Path(name).parts != ("networks", f"{row.Method}__{row.Case}.csv"):
            fail(f"invalid network identity: {row.Network}")
        if not (STAGE / name).is_file() or sha256(STAGE / name) != row.NetworkSHA256:
            fail(f"missing network or hash mismatch: {name}")
        hashes[name] = row.NetworkSHA256
    return runtime, hashes


def main():
    runtime, hashes = validate_batch()
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".beeline-promotion-", dir=TARGET.parent) as temporary:
        candidate = Path(temporary) / "results"
        candidate.mkdir()
        for name, digest in hashes.items():
            if name == "status.json":
                continue
            destination = candidate / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(STAGE / name, destination)
            if sha256(destination) != digest:
                fail(f"batch changed during copying: {name}")
        columns = ["Case", "Method", "DatasetType", "RuntimeSeconds", "TimingScope"]
        runtime[columns].to_csv(candidate / "runtime/runtime_all_methods.tsv", sep="\t", index=False)
        provenance = runtime[columns].assign(Origin=STAGE.name, Status="complete", Error="")
        provenance.to_csv(candidate / "runtime/runtime_provenance.tsv", sep="\t", index=False)
        backup = publish_directory(candidate, TARGET, BACKUP_ROOT, preserve=KEEP_IN_PLACE)
    print(f"promoted {len(runtime)} records from {STAGE} into {TARGET}")
    if backup:
        print(f"previous formal results verified in {backup}")


if __name__ == "__main__":
    main()
