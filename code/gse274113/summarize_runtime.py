from pathlib import Path
import csv
import tempfile

ROOT = Path(__file__).resolve().parents[2]
MEASUREMENTS = ROOT / "results/gse274113/runtime"


def read_table(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize_runtime(measurements=MEASUREMENTS, output=ROOT / "results/gse274113/runtime", comparators=MEASUREMENTS / "comparator_timings.tsv"):
    measurements, output = Path(measurements), Path(output)
    try:
        timing_source = measurements.resolve().relative_to(ROOT)
    except ValueError:
        timing_source = measurements.resolve()
    timings = []
    configurations = (("optimized_1core", 1, "single core"), ("optimized_8core", 8, "8 cores (fork)"))
    for configuration, _, _ in configurations:
        for day in (7, 9, 11, 14):
            rows = read_table(measurements / f"{configuration}_day_{day}_timing.tsv")
            timings.append(dict(configuration=configuration, **rows[0]))
    comparison = []
    for configuration, cores, label in configurations:
        selected = [row for row in timings if row["configuration"] == configuration]
        seconds = sum(float(row["inference_export_seconds"]) for row in selected)
        wrapper = sum(float(row["wrapper_seconds"]) for row in selected)
        comparison.append(dict(method="multiCSN", core_configuration=label, cores_measured=cores,
            inference_seconds=f"{seconds:.3f}", scope="timepoint normalization, candidate preparation, three-layer fitting, checkpoints, chain projection, export and validation; four timepoints",
            timing_source=f"{timing_source}/*_timing.tsv", cited_in="Figure S7a",
            note=f"Full wrapper including input reading and initial identity checks: {wrapper:.3f} s."))
    comparison.extend(read_table(comparators))
    tables = {"multiomic_runtime_comparison.tsv": comparison,
              "multiomic_runtime_timepoints.tsv": timings}
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".runtime-tables-", dir=output) as temporary:
        for name, rows in tables.items():
            path = Path(temporary) / name
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, rows[0], delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            path.replace(output / name)
    return tables


if __name__ == "__main__":
    summarize_runtime()
