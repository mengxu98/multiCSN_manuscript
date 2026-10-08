#!/usr/bin/env python3

import os
import time
from pathlib import Path


def process_snapshot():
    parents = {}
    rss = {}
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            stat = (entry / "stat").read_text()
            tail = stat[stat.rfind(")") + 2 :].split()
            parents[int(entry.name)] = int(tail[1])
            pages = int((entry / "statm").read_text().split()[1])
            rss[int(entry.name)] = pages * os.sysconf("SC_PAGE_SIZE") // 1024
        except (FileNotFoundError, PermissionError, ProcessLookupError, ValueError, IndexError):
            continue
    return parents, rss


def descendants(root, parents):
    selected = {root}
    changed = True
    while changed:
        changed = False
        for pid, parent in parents.items():
            if parent in selected and pid not in selected:
                selected.add(pid)
                changed = True
    return selected


def sample_rss(root_pid, output, interval=0.25):
    peak_rss = 0
    peak_processes = 0
    samples = 0
    while True:
        parents, rss = process_snapshot()
        selected = descendants(root_pid, parents)
        live = [pid for pid in selected if pid in rss]
        if live:
            current_rss = sum(rss[pid] for pid in live)
            if current_rss > peak_rss:
                peak_rss = current_rss
                peak_processes = len(live)
            samples += 1
        if root_pid not in parents and root_pid not in rss:
            break
        time.sleep(interval)

    Path(output).write_text(
        "process_tree_rss_sum_peak_kb\t{}\n"
        "process_tree_processes_at_peak\t{}\n"
        "process_tree_rss_samples\t{}\n"
        "process_tree_rss_interval_seconds\t{}\n".format(
            peak_rss, peak_processes, samples, interval
        )
    )
