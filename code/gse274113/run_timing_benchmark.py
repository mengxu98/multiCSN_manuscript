#!/usr/bin/env python3
"""Measure multimodal method runtime on the four GSE274113 control time points.

Protocol (documented in the emitted protocol.md):
  - one host (131), one method at a time, single inference core for the methods
    that expose a core setting;
  - identical per-time-point control inputs taken from the candidate catalogue;
  - every run is wrapped in `/usr/bin/time -v`, recording wall/user/system
    seconds, percent CPU and peak RSS;
  - outputs go to an isolated benchmark tree under runs/, never to results/.

The method and time-point sweep settings are declared in this script.
"""
from __future__ import annotations

import os
import re
import shlex
import subprocess
import time
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code/functions"))
from runtime import r_function_command, python_function_command, rscript_executable, PYTHON_LOCATIONS, COMPARATOR_R_LIBRARY, task_environment
OUTPUT_ROOT = REPO / "runs" / "gse274113" / "timing_benchmark"
LOG_ROOT = OUTPUT_ROOT / "logs"
TIMING_TSV = OUTPUT_ROOT / "timing.tsv"
RSCRIPT = rscript_executable()
LINGER_PYTHON = str(PYTHON_LOCATIONS["linger"][0])
LINGER_BIN = str(Path(LINGER_PYTHON).parent)
SCENICPLUS_PYTHON = str(PYTHON_LOCATIONS["scenicplus"][0])
CANDIDATES = REPO / "data/gse274113/objects/GSE274113_control_full/candidates"
TIME_POINTS = {
    "day 7": ("day_7", 3304),
    "day 9": ("day_9", 4167),
    "day 11": ("day_11", 1903),
    "day 14": ("day_14", 2073),
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


def run_step(method: str, timepoint: str, label: str, command: list[str], env_extra: dict) -> None:
    slug = TIME_POINTS[timepoint][0] if timepoint in TIME_POINTS else timepoint.replace(" ", "_")
    log_path = LOG_ROOT / f"{method}__{slug}.time.log"
    if log_path.exists():
        previous = log_path.read_text(errors="replace")
        if "Exit status: 0" in previous:
            print(f"SKIP completed {log_path.name}", flush=True)
            return
        failed = log_path.with_suffix(f".failed_{time.strftime('%Y%m%dT%H%M%S')}")
        log_path.rename(failed)
        print(f"RETRY {log_path.name} (previous attempt failed; kept {failed.name})", flush=True)
    load_start = os.getloadavg()
    env = task_environment()
    env.update(env_extra)
    started = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    with log_path.open("w") as log:
        completed = subprocess.run(
            ["/usr/bin/time", "-v", *command], stdout=log, stderr=subprocess.STDOUT,
            env=env, cwd=REPO,
        )
    text = log_path.read_text(errors="replace")
    wall = parse_elapsed(ELAPSED.search(text).group("value")) if ELAPSED.search(text) else float("nan")
    user = float(USER.search(text).group("value")) if USER.search(text) else float("nan")
    system = float(SYS.search(text).group("value")) if SYS.search(text) else float("nan")
    cpu = float(CPU.search(text).group("value")) if CPU.search(text) else float("nan")
    rss = int(RSS.search(text).group("value")) if RSS.search(text) else 0
    load_end = os.getloadavg()
    header = not TIMING_TSV.exists()
    with TIMING_TSV.open("a") as stream:
        if header:
            stream.write(
                "method\ttimepoint\tstage\twall_seconds\tuser_seconds\tsystem_seconds\t"
                "percent_cpu\tpeak_rss_kb\tload_1min_start\tload_1min_end\tstarted\t"
                "exit_status\tcommand\n"
            )
        stream.write(
            f"{method}\t{timepoint}\t{label}\t{wall:.3f}\t{user:.2f}\t{system:.2f}\t{cpu:.0f}\t"
            f"{rss}\t{load_start[0]:.2f}\t{load_end[0]:.2f}\t{started}\t{completed.returncode}\t"
            f"{shlex.join(command)}\n"
        )
    print(f"{method} {timepoint} {label}: wall={wall:.1f}s rss={rss/1048576:.1f}GiB exit={completed.returncode}",
          flush=True)
    if completed.returncode != 0:
        raise SystemExit(f"{method} {timepoint} failed; inspect {log_path}")


def multi_csn_steps(timepoints):
    for timepoint in timepoints:
        slug, cells = TIME_POINTS[timepoint]
        output = OUTPUT_ROOT / "multiCSN" / slug
        if output.exists():
            print(f"SKIP existing {output}", flush=True)
            continue
        yield timepoint, "inference_export_multilayer", r_function_command(
            "code/functions/network_fit.R", "fit_state_networks", executable=RSCRIPT,
            dataset="GSE274113", input=CANDIDATES / "GSE274113_NT_timepoint_candidates.rds",
            validation=CANDIDATES / "candidate_validation.tsv", output=output,
            state=timepoint, expected_cells=cells, response_chunk=1024, profile="gse274113"), {}



def pando_steps(timepoints):
    output = OUTPUT_ROOT / "Pando"
    if output.exists():
        print(f"SKIP existing {output}", flush=True)
        return
    yield "all timepoints", "inference_export_all_timepoints", r_function_command(
        "code/gse274113/run_pando.R", "run_pando", executable=RSCRIPT, output_root=output), {
        "R_LIBS_USER": str(COMPARATOR_R_LIBRARY),
    }



def linger_steps(timepoints):
    for timepoint in timepoints:
        slug = TIME_POINTS[timepoint][0]
        output = OUTPUT_ROOT / "LINGER"
        marker = output / slug / "networks"
        if marker.exists():
            print(f"SKIP existing {marker}", flush=True)
            continue
        yield timepoint, "prepare_preprocess_train_networks_export", python_function_command(
            "code/gse274113/run_linger.py", "run_timepoint", executable=LINGER_PYTHON,
            timepoint=slug, output_root=output), {
            "PYTHONHASHSEED": "42", "PATH": LINGER_BIN + os.pathsep + os.environ.get("PATH", ""),
        }



def scenicplus_steps(timepoints):
    # Scope A uses the archived production pycisTopic models and needs no
    # auxiliary resources; scope B (prepare + LDA) additionally requires
    # SCENICplus_auxiliary/{hg38-blacklist.v2.bed,Mallet-202108}, which are
    # currently missing from data/ and from every backup.
    production_object = REPO / "data/gse274113/objects/GSE274113_SCENICplus_full"
    if not (production_object / "cistopic_models").exists():
        raise SystemExit("archived pycisTopic models are missing; cannot time scope A")
    yield "all timepoints", "scope_A_downstream_eregulons", python_function_command(
        "code/gse274113/run_scenicplus.py", "run_timing", executable=SCENICPLUS_PYTHON), {}
    print("scope B (prepare + pycisTopic LDA) not started: "
          "SCENICplus_auxiliary resources are missing", flush=True)


def main() -> None:
    methods = METHODS
    timepoints = tuple(TIME_POINTS)
    for method in methods:
        if method not in METHODS:
            raise SystemExit(f"Unknown method: {method}")
    for timepoint in timepoints:
        if timepoint not in TIME_POINTS:
            raise SystemExit(f"Unknown timepoint: {timepoint}")
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    LOG_ROOT.mkdir(parents=True, exist_ok=True)
    protocol = OUTPUT_ROOT / "protocol.md"
    if not protocol.exists():
        protocol.write_text(
            "# GSE274113 multimodal method timing benchmark\n\n"
            "- Host: 131 (Intel Xeon E5-2680 v3), methods run one at a time.\n"
            "- Inputs: the same per-time-point control candidate catalogue for every method.\n"
            "- Scope: network inference and export. SCENIC+ starts from its\n"
            "  archived pycisTopic models (scope A) and additionally records the full\n"
            "  LDA-inclusive pipeline (scope B); multiCSN/Pando/LINGER run their method-specific\n"
            "  pipelines because those have no separately archived upstream model.\n"
            "- Single inference core for multiCSN/Pando (the LINGER workflow is also\n"
            "  single-process); SCENIC+ uses its configured worker count and the row\n"
            "  records that value.\n"
            "- Each run is wrapped in `/usr/bin/time -v`; timing.tsv records wall,\n"
            "  user and system seconds, percent CPU, peak RSS and host load.\n"
        )
    generators = {
        "multiCSN": multi_csn_steps,
        "Pando": pando_steps,
        "LINGER": linger_steps,
        "SCENICplus": scenicplus_steps,
    }
    for method in methods:
        generator = generators[method]
        for timepoint, label, command, env_extra in generator(timepoints):
            run_step(method, timepoint, label, command, env_extra)
    print("benchmark batch finished", flush=True)


if __name__ == "__main__":
    main()
