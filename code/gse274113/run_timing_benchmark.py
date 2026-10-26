#!/usr/bin/env python3
"""Measure multimodal method runtime on the four GSE274113 control time points.

  - one host, one method at a time, single inference core for the methods
    that expose a core setting;
  - identical per-time-point control inputs taken from the candidate catalogue;
  - every run is wrapped in `/usr/bin/time -v`, recording wall/user/system
    seconds, percent CPU and peak RSS;
  - outputs go to an isolated benchmark tree under runs/, never to results/.

The method and time-point sweep settings are declared in this script.
"""
from __future__ import annotations

import os
import math
import re
import csv
import tempfile
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code/functions"))
from runtime import r_function_command, python_function_command, rscript_executable, PYTHON_LOCATIONS, task_environment
OUTPUT_ROOT = REPO / "runs" / "gse274113" / "timing_benchmark"
TIMING_TSV = OUTPUT_ROOT / "timing.tsv"
RSCRIPT = rscript_executable()
LINGER_PYTHON = str(PYTHON_LOCATIONS["linger"])
LINGER_BIN = str(Path(LINGER_PYTHON).parent)
SCENICPLUS_PYTHON = str(PYTHON_LOCATIONS["scenicplus"])
CANDIDATES = REPO / "data/gse274113/objects/GSE274113_control_full/candidates"
TIME_POINTS = {
    "day 7": "day_7",
    "day 9": "day_9",
    "day 11": "day_11",
    "day 14": "day_14",
}
METHODS = ("multiCSN", "Pando", "LINGER", "SCENICplus")
ELAPSED = re.compile(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (?P<value>[\d:.]+)")
USER = re.compile(r"User time \(seconds\): (?P<value>[\d.]+)")
SYS = re.compile(r"System time \(seconds\): (?P<value>[\d.]+)")
CPU = re.compile(r"Percent of CPU this job got: (?P<value>\d+)%")
RSS = re.compile(r"Maximum resident set size \(kbytes\): (?P<value>\d+)")


def parse_elapsed(text: str) -> float:
    parts = text.split(":")
    seconds = float(parts[-1])
    if len(parts) >= 2:
        seconds += 60 * float(parts[-2])
    if len(parts) >= 3:
        seconds += 3600 * float(parts[-3])
    return seconds


def timing_complete(row, output, marker):
    try:
        wall = float(row["wall_seconds"])
        if row["exit_status"] != "0" or not math.isfinite(wall) or wall < 0:
            return False
        if not (output / marker).is_file():
            return False
        if row["method"] == "LINGER":
            with (output / "phase_timings.tsv").open() as handle:
                phases = list(csv.DictReader(handle, delimiter="\t"))
            expected = {"prepare", "preprocess", "train", "networks", "export"}
            return len(phases) == len(expected) and {r["phase"] for r in phases} == expected and all(
                r["timepoint"].replace(" ", "_") == row["timepoint"].replace(" ", "_")
                and math.isfinite(float(r["seconds"])) and float(r["seconds"]) >= 0 for r in phases)
        return True
    except (OSError, KeyError, ValueError, TypeError):
        return False


def run_step(method: str, timepoint: str, label: str, command: list[str], env_extra: dict,
             output: Path, marker: str) -> None:
    previous = None
    if TIMING_TSV.exists():
        with TIMING_TSV.open() as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                if row["method"] == method and row["timepoint"] == timepoint and row["stage"] == label:
                    previous = row
    if previous is not None and timing_complete(previous, output, marker):
        print(f"SKIP completed {method} {timepoint} {label}", flush=True)
        return
    if output.exists() or (previous is not None and previous["exit_status"] == "0"):
        raise RuntimeError(f"Incomplete timing record or outputs for {method} {timepoint}: {output}. "
                           "Keep the existing results and use a new timing output directory to remeasure.")
    env = task_environment()
    env.update(env_extra)
    with tempfile.TemporaryDirectory(prefix="multicsn-timing-") as temporary:
        timing = Path(temporary) / "time.txt"
        completed = subprocess.run(
            ["/usr/bin/time", "-v", "-o", str(timing), *command], env=env, cwd=REPO,
        )
        text = timing.read_text(errors="replace")
    wall = parse_elapsed(ELAPSED.search(text).group("value")) if ELAPSED.search(text) else float("nan")
    user = float(USER.search(text).group("value")) if USER.search(text) else float("nan")
    system = float(SYS.search(text).group("value")) if SYS.search(text) else float("nan")
    cpu = float(CPU.search(text).group("value")) if CPU.search(text) else float("nan")
    rss = int(RSS.search(text).group("value")) if RSS.search(text) else 0
    record = dict(method=method, timepoint=timepoint, stage=label,
                  wall_seconds=f"{wall:.3f}", user_seconds=f"{user:.2f}", system_seconds=f"{system:.2f}",
                  percent_cpu=f"{cpu:.0f}", peak_rss_kb=rss, exit_status=completed.returncode)
    header = not TIMING_TSV.exists()
    if not header:
        with TIMING_TSV.open() as stream:
            fields = next(csv.reader(stream, delimiter="\t"))
    else:
        fields = list(record)
    with TIMING_TSV.open("a", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t")
        if header:
            writer.writeheader()
        writer.writerow(record)
    print(f"{method} {timepoint} {label}: wall={wall:.1f}s rss={rss/1048576:.1f}GiB exit={completed.returncode}",
          flush=True)
    if completed.returncode != 0:
        raise SystemExit(f"{method} {timepoint} failed with exit status {completed.returncode}")
    if not timing_complete({key: str(value) for key, value in record.items()}, output, marker):
        raise RuntimeError(f"Missing valid timing values or final outputs for {method} {timepoint}")


def multi_csn_steps(timepoints):
    for timepoint in timepoints:
        slug = TIME_POINTS[timepoint]
        output = OUTPUT_ROOT / "multiCSN" / slug
        yield timepoint, "inference_export_multilayer", r_function_command(
            "code/functions/network_fit.R", "fit_state_networks", executable=RSCRIPT,
            dataset="GSE274113", input=CANDIDATES / "GSE274113_NT_timepoint_candidates.rds",
            output=output,
            state=timepoint, response_chunk=1024, profile="gse274113"), {}, output, "endpoint_summary.tsv"



def pando_steps(timepoints):
    for cores in (1, 8):
        output = OUTPUT_ROOT / ("Pando" if cores == 1 else "Pando_8core")
        timepoint = "all timepoints" if cores == 1 else "all timepoints 8core"
        yield timepoint, f"inference_export_{cores}core_all_timepoints", r_function_command(
            "code/gse274113/run_pando.R", "run_pando", executable=RSCRIPT, cores=cores, output_root=output), {}, output, "timepoint_run_summary.tsv"



def linger_steps(timepoints):
    for timepoint in timepoints:
        slug = TIME_POINTS[timepoint]
        output = OUTPUT_ROOT / "LINGER"
        yield timepoint, "prepare_preprocess_train_networks_export", python_function_command(
            "code/gse274113/run_linger.py", "run_timepoint", executable=LINGER_PYTHON,
            timepoint=slug, output_root=output), {
            "PYTHONHASHSEED": "42", "PATH": LINGER_BIN + os.pathsep + os.environ.get("PATH", ""),
        }, output / slug, "networks/network_summary.tsv"



def scenicplus_steps(timepoints):
    yield "all timepoints", "scope_A_downstream_eregulons", python_function_command(
        "code/gse274113/run_scenicplus.py", "run_timing", executable=SCENICPLUS_PYTHON), {}, (
            OUTPUT_ROOT / "SCENICplus/output/SCENICplus"), "timepoint_run_summary.tsv"
    print("scope B (prepare + pycisTopic LDA) not started: "
          "SCENICplus_auxiliary resources are missing", flush=True)


def main() -> None:
    methods = METHODS
    timepoints = tuple(TIME_POINTS)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    generators = {
        "multiCSN": multi_csn_steps,
        "Pando": pando_steps,
        "LINGER": linger_steps,
        "SCENICplus": scenicplus_steps,
    }
    for method in methods:
        generator = generators[method]
        for timepoint, label, command, env_extra, output, marker in generator(timepoints):
            run_step(method, timepoint, label, command, env_extra, output, marker)
    print("benchmark batch finished", flush=True)


if __name__ == "__main__":
    main()
