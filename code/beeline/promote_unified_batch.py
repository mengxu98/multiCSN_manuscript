#!/usr/bin/env python3
import math
import shutil
import tempfile
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "functions"))
from runtime import use_python
from file_io import publish_directory
sys.path.insert(0, str(Path(__file__).absolute().parent))
from run_unified_benchmark import METHODS, case_paths
if __name__ == "__main__":
    use_python("beeline")
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
STAGE = REPO / "runs/beeline/unified_formal"
TARGET = REPO / "results/beeline"
KEEP_IN_PLACE = {"core_scaling"}
EXPECTED_METHODS = len(METHODS)
EXPECTED_CASES = 500
EXPECTED_RECORDS = EXPECTED_METHODS * EXPECTED_CASES
METRIC_FILES = (
    "curated__primary_metrics.csv", "curated__other_metrics.csv",
    "synthetic__primary_metrics.csv", "synthetic__other_metrics.csv",
    "combined__other_metrics.csv",
    "case_metrics.tsv", "case_ranks.tsv", "rank_summary.tsv",
)


def validate_batch():
    runtime = pd.read_csv(STAGE / "runtime/case_runtime.tsv", sep="\t")
    required = {"Case", "Method", "DatasetType", "RuntimeSeconds", "TimingScope"}
    if not required.issubset(runtime.columns) or len(runtime) != EXPECTED_RECORDS:
        raise ValueError("Incomplete benchmark runtimes")
    if runtime.duplicated(["Case", "Method"]).any() or set(runtime.Method) != set(METHODS):
        raise ValueError("Duplicate or unexpected benchmark cases")
    cases = {path.parent.name for path in case_paths(0)}
    if len(cases) != EXPECTED_CASES or any(set(group.Case) != cases for _, group in runtime.groupby("Method")):
        raise ValueError("Incomplete case coverage")
    if not runtime.RuntimeSeconds.map(lambda value: math.isfinite(value) and value > 0).all():
        raise ValueError("Invalid benchmark runtimes")
    if not runtime.TimingScope.eq("inference_and_export").all():
        raise ValueError("Inconsistent timing scope")
    files = [f"networks/{row.Method}__{row.Case}.csv" for row in runtime.itertuples()]
    files += [f"metrics/{name}" for name in METRIC_FILES]
    files += ["runtime/case_runtime.tsv", "runtime/runtime_summary.tsv"]
    files += [str(path.relative_to(STAGE)) for path in (STAGE / "evaluation").glob("*.tsv")]
    if not (STAGE / "evaluation/case_metrics.tsv").is_file():
        raise FileNotFoundError("Missing benchmark evaluation")
    for name in files:
        if not (STAGE / name).is_file():
            raise FileNotFoundError(name)
    return runtime, files


def main():
    runtime, files = validate_batch()
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".beeline-", dir=TARGET.parent) as temporary:
        candidate = Path(temporary) / "results"
        candidate.mkdir()
        for name in files:
            destination = candidate / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(STAGE / name, destination)
        runtime[["Case", "Method", "DatasetType", "RuntimeSeconds", "TimingScope"]].to_csv(
            candidate / "runtime/runtime_all_methods.tsv", sep="\t", index=False)
        publish_directory(candidate, TARGET, preserve=KEEP_IN_PLACE)
    shutil.rmtree(STAGE)
    print(f"Saved {len(runtime)} benchmark records in {TARGET}")


if __name__ == "__main__":
    main()
