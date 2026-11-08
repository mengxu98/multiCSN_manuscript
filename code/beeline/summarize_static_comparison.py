import csv
import hashlib
import json
import platform
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import statistics
import sys

repo = Path(__file__).resolve().parents[2]
run = Path(sys.argv[1]).resolve()
destination = Path(sys.argv[2]).resolve()
preflight = Path(sys.argv[3]).resolve()
if not (run / "COMPLETE.txt").is_file():
    raise RuntimeError("The paired experiment is incomplete")
if destination.exists():
    raise FileExistsError(destination)


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(name, rows):
    with (destination / name).open("x") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


frozen = json.loads(preflight.read_text())
for item in frozen["files"]:
    path = repo / item["path"]
    if path.stat().st_size != item["bytes"] or hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
        raise RuntimeError(f"Frozen input or baseline changed: {path}")
replay = read_table(run / "baseline_replay_checks.tsv")
if len(replay) != 500 or any(row["NetworkMatch"] != "TRUE" or
                           float(row["ReplayMetricMaxError"]) >= 1e-10 or
                           float(row["BaselineMetricMaxError"]) >= 1e-10 for row in replay):
    raise RuntimeError("Historical replay verification failed")
rows = read_table(run / "all_method_case_metrics.tsv")
if len(rows) != 6000 or len({(row["Case"], row["Method"]) for row in rows}) != 6000:
    raise RuntimeError("Expected 6000 unique case/method entries")
destination.mkdir(parents=True)
for name in ["cases.tsv", "mode_case_metrics.tsv", "baseline_replay_checks.tsv",
             "all_method_case_metrics.tsv", "method_summary.tsv", "group_summary.tsv",
             "paired_differences.tsv", "paired_summary.tsv", "case_runtime.tsv",
             "runtime_repeats.tsv", "runtime_summary.tsv", "sessionInfo.txt", "COMPLETE.txt"]:
    (destination / name).write_bytes((run / name).read_bytes())
metrics = ["AP", "AUPRC", "AUROC", "EarlyPrecision", "EPR", "ACC", "Precision", "Recall", "F1"]
groups = defaultdict(list)
for row in rows:
    prefix = row["DatasetGroup"] + "-"
    if not row["Case"].startswith(prefix):
        raise RuntimeError(f"Case/model mismatch: {row['Case']}")
    fields = row["Case"][len(prefix):].split("-")
    if len(fields) not in (2, 3) or any(not value.isdigit() for value in fields):
        raise RuntimeError(f"Unrecognized configuration: {row['Case']}")
    cells, replicate = fields[:2]
    dropout = fields[2] if len(fields) == 3 else "0"
    if not 1 <= int(replicate) <= 10:
        raise RuntimeError(f"Unexpected replicate: {row['Case']}")
    groups[(row["DatasetType"], row["DatasetGroup"], cells,
            dropout, row["Method"])].append(row)
configurations = []
for key, values in sorted(groups.items()):
    if len(values) != 10:
        raise RuntimeError(f"Configuration is not complete: {key}")
    item = dict(zip(["DatasetType", "Model", "Cells", "DropoutPercent", "Method"], key))
    item["Cases"] = len(values)
    item.update({metric: statistics.mean(float(row[metric]) for row in values) for metric in metrics})
    configurations.append(item)
write_table("configuration_summary.tsv", configurations)
scenarios = []
for mode in ["static", "pseudotime"]:
    selected = [row for row in rows if not row["Method"].startswith("multiCSN") or
                row["Method"] == f"multiCSN ({mode})"]
    ranks = {}
    by_case = defaultdict(list)
    for row in selected:
        by_case[row["Case"]].append(row)
    for case, values in by_case.items():
        if len(values) != 11:
            raise RuntimeError(f"Expected eleven entries in scenario {mode}, {case}")
        for metric in metrics:
            ordered = sorted(values, key=lambda row: -float(row[metric]))
            for row in ordered:
                value = float(row[metric])
                positions = [i + 1 for i, other in enumerate(ordered) if float(other[metric]) == value]
                ranks[(case, row["Method"], metric)] = statistics.mean(positions)
    for scope in ["all", "curated", "synthetic"]:
        grouped = defaultdict(list)
        for row in selected:
            if scope == "all" or row["DatasetType"] == scope:
                grouped[row["Method"]].append(row)
        for method, values in sorted(grouped.items()):
            item = {"Scenario": mode, "Scope": scope, "Method": method, "Cases": len(values)}
            for metric in metrics:
                item[metric] = statistics.mean(float(row[metric]) for row in values)
                item[f"{metric}_rank11"] = statistics.mean(ranks[(row["Case"], method, metric)] for row in values)
            scenarios.append(item)
write_table("scenario_summary.tsv", scenarios)
provenance = dict(frozen)
provenance["host"] = {"system": platform.system(), "architecture": platform.machine(),
                      "macOS": platform.mac_ver()[0]}
provenance["completed_utc"] = datetime.now(timezone.utc).isoformat()
provenance["validation"] = {
    "cases": 500, "baseline_network_matches": 500, "timed_calls": 5000,
    "joint_method_entries": 12, "joint_case_rows": 6000,
    "configuration_method_rows": len(configurations), "frozen_files_unchanged": len(frozen["files"]),
    "maximum_baseline_metric_error": max(float(row["BaselineMetricMaxError"]) for row in replay),
    "runtime_comparison": "Five alternating same-host repetitions per mode; historical comparator timings are not rerun or ranked with these timings."
}
provenance["outputs"] = [{"path": str(path.relative_to(repo)), "bytes": path.stat().st_size,
                          "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                         for path in sorted(destination.iterdir())]
provenance["summary_code_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
package_archive = preflight.parent / "inferCSN_installed.tar.gz"
if package_archive.is_file():
    provenance["installed_package_archive"] = {
        "path": str(package_archive.relative_to(repo)), "bytes": package_archive.stat().st_size,
        "sha256": hashlib.sha256(package_archive.read_bytes()).hexdigest()
    }
elif "installed_package_archive" in frozen:
    provenance["installed_package_archive"] = dict(frozen["installed_package_archive"],
                                                   availability="recorded_hash_only")
else:
    provenance["installed_package_archive"] = {"availability": "not_supplied"}
provenance["configuration_labels"] = "Derived from frozen case and model identifiers; no external case manifest is used."
(destination / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
print(json.dumps(provenance["validation"], indent=2))
