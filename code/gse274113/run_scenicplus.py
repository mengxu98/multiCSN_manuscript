#!/usr/bin/env python3
"""Run four independent SCENIC+ timepoint networks with validated offline inputs."""
from __future__ import annotations
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "functions"))
from runtime import (REPO, EXECUTION_ROOT, python_executable, use_python,
                     run_task, task_environment, tee_log, utc_now, sha256)

CORES = 12
TIMEPOINTS = [("day_7", "day 7", 3304), ("day_9", "day 9", 4167),
              ("day_11", "day 11", 1903), ("day_14", "day 14", 2073)]
OBJECT_ROOT = REPO / "data/gse274113/objects/GSE274113_SCENICplus_full"
OUTPUT_PARENT = REPO / "runs/gse274113/method_fits"
LOG_ROOT = REPO / "logs/gse274113/comparators"
OFFLINE_PRODUCER = REPO / "code/gse274113/prepare_scenicplus_offline_genome_resources.py"
TIMING_ROOT = REPO / "runs/gse274113/timing_benchmark/SCENICplus"
TASK_SETTINGS = {
    "network": dict(object_root=OBJECT_ROOT, output_parent=OUTPUT_PARENT,
                    log_root=LOG_ROOT, offline_producer=OFFLINE_PRODUCER,
                    cores=CORES, monitor=True),
    "timing": dict(object_root=OBJECT_ROOT, output_parent=TIMING_ROOT / "output",
                   log_root=TIMING_ROOT / "logs",
                   offline_producer=TIMING_ROOT / "provenance/prepare_scenicplus_offline_genome_resources.py",
                   cores=CORES, monitor=False),
}

def validate_offline_genome_resources(settings):
    import csv
    from pathlib import Path

    root = Path(settings['offline_resource_root']).resolve()
    producer = Path(settings['offline_resource_producer']).resolve()


    with (root / "artifact_manifest.tsv").open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    if not manifest or list(manifest[0]) != ["file", "bytes", "sha256"]:
        raise RuntimeError("Invalid offline genome artifact manifest schema")
    expected_files = {
        "genome_annotation.tsv", "chromsizes.tsv", "source_manifest.tsv",
        "resource_validation.tsv",
    }
    if {row["file"] for row in manifest} != expected_files:
        raise RuntimeError("Offline genome artifact manifest has incomplete coverage")
    for row in manifest:
        path = root / row["file"]
        if not path.is_file() or path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"]:
            raise RuntimeError(f"Stale offline genome artifact: {path}")

    with (root / "resource_validation.tsv").open(newline="") as handle:
        validation = {row["field"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
    required = {
        "assembly": "GRCh38.p14",
        "annotation_release": "Ensembl 116",
        "chromosome_style": "UCSC",
        "primary_chromosomes": "25",
        "producer_script_sha256": sha256(producer),
    }
    for field, value in required.items():
        if validation.get(field) != value:
            raise RuntimeError(f"Offline genome validation mismatch for {field}")

    with (root / "genome_annotation.tsv").open(newline="") as handle:
        annotation = csv.DictReader(handle, delimiter="\t")
        expected = [
            "Chromosome", "Start", "End", "Strand", "Gene",
            "Transcription_Start_Site", "Transcript_type",
        ]
        if annotation.fieldnames != expected or sum(1 for _ in annotation) != int(validation["annotation_rows"]):
            raise RuntimeError("Offline genome annotation schema or row count mismatch")
    with (root / "chromsizes.tsv").open(newline="") as handle:
        chromsizes = list(csv.DictReader(handle, delimiter="\t"))
    if list(chromsizes[0]) != ["Chromosome", "Start", "End"] or len(chromsizes) != 25:
        raise RuntimeError("Offline chromsizes schema or row count mismatch")
    expected_chromosomes = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
    if {row["Chromosome"] for row in chromsizes} != expected_chromosomes:
        raise RuntimeError("Offline chromsizes does not cover the 25 primary chromosomes")
    if any(int(row["Start"]) != 0 or int(row["End"]) <= 0 for row in chromsizes):
        raise RuntimeError("Offline chromsizes contains invalid bounds")
    print("offline_genome_resource_validation\tPASS")


def validate_timepoint_artifacts(settings):
    from pathlib import Path

    import pandas as pd

    root = Path(settings['time_work']).resolve()
    manifest_path = root / "artifact_manifest.tsv"
    input_manifest_path = root / "input_manifest.tsv"
    validation_path = root / "run_validation.tsv"
    for path in (manifest_path, input_manifest_path, validation_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing non-empty resumable SCENIC+ record: {path}")


    def validate_manifest(path, base, require_relative):
        table = pd.read_csv(path, sep="\t", dtype=str)
        if list(table.columns) != ["file", "bytes", "sha256"] or table.empty or table["file"].duplicated().any():
            raise RuntimeError(f"Invalid SCENIC+ hash manifest schema: {path}")
        resolved = []
        for value in table["file"]:
            candidate = Path(value)
            if require_relative:
                if candidate.is_absolute() or ".." in candidate.parts:
                    raise RuntimeError(f"Unsafe relative artifact path in {path}: {value}")
                candidate = (base / candidate).resolve()
                if candidate.parent != base and base not in candidate.parents:
                    raise RuntimeError(f"Artifact escaped its timepoint directory: {value}")
            else:
                if not candidate.is_absolute():
                    raise RuntimeError(f"Input manifest path is not absolute: {value}")
                candidate = candidate.resolve()
            resolved.append(candidate)
        for row, candidate in zip(table.itertuples(index=False), resolved):
            if not candidate.is_file() or candidate.stat().st_size != int(row.bytes) or sha256(candidate) != row.sha256:
                raise RuntimeError(f"Stale or corrupt SCENIC+ manifest member: {candidate}")
        return table, resolved

    artifact, _ = validate_manifest(manifest_path, root, True)
    actual = sorted(
        str(path.relative_to(root))
        for path in root.rglob("*")
        if path.is_file() and path != manifest_path and "tmp" not in path.relative_to(root).parts
    )
    if set(actual) != set(artifact["file"]):
        raise RuntimeError("SCENIC+ timepoint manifest does not exactly cover the resumable directory")
    inputs, _ = validate_manifest(input_manifest_path, root, False)
    region_root = Path(settings['region_sets']).resolve()
    region_files = sorted(path.resolve() for path in region_root.rglob("*") if path.is_file())
    if not region_files:
        raise RuntimeError("SCENIC+ region_sets directory is empty")
    expected_inputs = {
        Path(settings[name]).resolve()
        for name in (
            "input_cistopic", "input_gex", "cfg_path", "resource_manifest",
            "snakefile", "default_config", "source_commit_file", "producer_script",
            "offline_genome_annotation", "offline_chromsizes", "offline_genome_manifest",
        )
    }
    expected_inputs.update(region_files)
    manifested_inputs = {Path(value).resolve() for value in inputs["file"]}
    if manifested_inputs != expected_inputs:
        raise RuntimeError("SCENIC+ input manifest does not exactly cover current model inputs")

    validation = pd.read_csv(validation_path, sep="\t", dtype=str)
    if list(validation.columns) != ["field", "value"] or validation["field"].duplicated().any():
        raise RuntimeError("SCENIC+ timepoint validation schema is invalid")
    values = validation.set_index("field")["value"]
    expected = {
        "timepoint": settings['timepoint'],
        "training_cells": str(int(settings['expected_cells'])),
        "network_pooling": "none",
        "cell_cap": "none",
        "gene_cap": "none",
        "region_cap": "none",
        "scenicplus_cores": str(int(settings['scenicplus_cores'])),
        "offline_genome_manifest_sha256": sha256(Path(settings['offline_genome_manifest']).resolve()),
        "producer_script_sha256": sha256(Path(settings['producer_script']).resolve()),
        "SCENICplus_source_commit": Path(settings['source_commit_file']).read_text().strip(),
    }
    for field, value in expected.items():
        if field not in values or values[field] != value:
            raise RuntimeError(f"SCENIC+ resumable validation mismatch for {field}")
    print(f"resume_artifact_validation\tPASS\t{settings['timepoint']}")


def write_config(settings):
    import yaml
    with open(settings['cfg_default']) as handle:
        config = yaml.safe_load(handle)
    config["input_data"].update({
        "cisTopic_obj_fname": settings['input_cistopic'],
        "GEX_anndata_fname": settings['input_gex'],
        "region_set_folder": settings['region_sets'],
        "ctx_db_fname": settings['ranking_db'],
        "dem_db_fname": settings['scores_db'],
        "path_to_motif_annotations": settings['motif_annotations'],
    })
    config["params_general"].update({
        "temp_dir": settings['work_tmp'],
        "n_cpu": int(settings['scenicplus_cores']),
        "seed": int(settings['seed']),
    })
    config["params_data_preparation"].update({
        "bc_transform_func": '"lambda x: f\'{x}\'"',
        "is_multiome": True,
        "key_to_group_by": "Timepoint",
    })
    with open(settings['cfg_path'], "w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)


def write_timepoint_records(settings):
    from pathlib import Path
    import pandas as pd

    root = Path(settings['time_work'])
    required = [
        "ACC_GEX.h5mu", "dem_results.hdf5", "ctx_results.hdf5",
        "cistromes_direct.h5ad", "cistromes_extended.h5ad", "tf_names.txt",
        "genome_annotation.tsv", "chromsizes.tsv", "search_space.tsv",
        "tf_to_gene_adj.tsv", "region_to_gene_adj.tsv", "eRegulon_direct.tsv",
        "eRegulons_extended.tsv", "AUCell_direct.h5mu", "AUCell_extended.h5mu",
        "scplusmdata.h5mu",
    ]
    missing = [name for name in required if not (root / name).is_file() or (root / name).stat().st_size == 0]
    if missing:
        raise RuntimeError("Missing SCENIC+ outputs: " + ", ".join(missing))

    def count_rows(path):
        with path.open("rt", errors="replace") as handle:
            return max(sum(1 for _ in handle) - 1, 0)


    counts = {name: count_rows(root / name) for name in [
        "tf_to_gene_adj.tsv", "region_to_gene_adj.tsv", "eRegulon_direct.tsv", "eRegulons_extended.tsv"
    ]}
    if any(value <= 0 for value in counts.values()):
        raise RuntimeError(f"One or more SCENIC+ edge tables are empty: {counts}")
    source_commit = Path(settings['source_commit_file']).read_text().strip()
    validation = pd.DataFrame({
        "field": [
            "dataset", "timepoint", "training_cells", "heldout_cells", "network_unit",
            "network_pooling", "tf_to_gene_rows", "region_to_gene_rows",
            "direct_eRegulon_rows", "extended_eRegulon_rows", "native_score_scope",
            "scenicplus_fit_seconds", "cell_cap", "gene_cap", "region_cap",
            "scenicplus_cores", "offline_genome_manifest_sha256",
            "method_native_region_set_selection",
            "SCENICplus_source_commit", "producer_script_sha256",
        ],
        "value": [
            "GSE274113", settings['timepoint'], int(settings['expected_cells']),
            126157, "independent control-timepoint eRegulon", "none",
            counts["tf_to_gene_adj.tsv"], counts["region_to_gene_adj.tsv"],
            counts["eRegulon_direct.tsv"], counts["eRegulons_extended.tsv"],
            "within_timepoint_native_endpoint_only", int(settings['fit_seconds']),
            "none", "none", "none",
            int(settings['scenicplus_cores']),
            sha256(Path(settings['offline_genome_manifest'])),
            "Otsu and top-3000 regions per selected topic; both collections retained",
            source_commit, sha256(Path(settings['producer_script'])),
        ],
    })
    validation.to_csv(root / "run_validation.tsv", sep="\t", index=False)

    input_paths = [
        Path(settings['input_cistopic']), Path(settings['input_gex']),
        Path(settings['cfg_path']), Path(settings['resource_manifest']),
        Path(settings['snakefile']), Path(settings['default_config']),
        Path(settings['source_commit_file']), Path(settings['producer_script']),
        Path(settings['offline_genome_annotation']),
        Path(settings['offline_chromsizes']),
        Path(settings['offline_genome_manifest']),
    ]
    region_set_paths = sorted(
        path for path in Path(settings['region_sets']).rglob("*") if path.is_file()
    )
    if not region_set_paths:
        raise RuntimeError("SCENIC+ region_sets directory contains no files")
    input_paths = list(dict.fromkeys(path.resolve() for path in input_paths + region_set_paths))
    pd.DataFrame({
        "file": [str(path) for path in input_paths],
        "bytes": [path.stat().st_size for path in input_paths],
        "sha256": [sha256(path) for path in input_paths],
    }).to_csv(root / "input_manifest.tsv", sep="\t", index=False)
    rows = []
    for path in sorted(p for p in root.rglob("*") if p.is_file() and "tmp" not in p.parts):
        if path.name == "artifact_manifest.tsv":
            continue
        rows.append({"file": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(rows).to_csv(root / "artifact_manifest.tsv", sep="\t", index=False)
    print(validation.to_string(index=False))


def write_suite_records(settings):
    from pathlib import Path
    import pandas as pd

    root = Path(settings['work_dir'])
    timepoints = [("day_7", "day 7", 3304), ("day_9", "day 9", 4167),
                  ("day_11", "day 11", 1903), ("day_14", "day 14", 2073)]
    rows = []
    for slug, timepoint, expected in timepoints:
        table = pd.read_csv(root / slug / "run_validation.tsv", sep="\t", dtype=str).set_index("field")["value"]
        if table["timepoint"] != timepoint or int(table["training_cells"]) != expected or table["network_pooling"] != "none":
            raise RuntimeError(f"SCENIC+ timepoint validation failed for {timepoint}")
        rows.append({
            "timepoint": timepoint,
            "training_cells": int(table["training_cells"]),
            "tf_to_gene_rows": int(table["tf_to_gene_rows"]),
            "region_to_gene_rows": int(table["region_to_gene_rows"]),
            "direct_eRegulon_rows": int(table["direct_eRegulon_rows"]),
            "extended_eRegulon_rows": int(table["extended_eRegulon_rows"]),
            "scenicplus_fit_seconds": float(table["scenicplus_fit_seconds"]),
            "scenicplus_cores": int(table["scenicplus_cores"]),
            "offline_genome_manifest_sha256": table["offline_genome_manifest_sha256"],
            "network_pooling": table["network_pooling"],
            "cell_cap": table["cell_cap"],
            "gene_cap": table["gene_cap"],
            "region_cap": table["region_cap"],
            "method_native_region_set_selection": table["method_native_region_set_selection"],
        })
    summary = pd.DataFrame(rows)
    if int(summary["training_cells"].sum()) != 11447 or len(summary) != 4:
        raise RuntimeError("SCENIC+ independent timepoint suite failed full-cell accounting.")
    summary.to_csv(root / "timepoint_run_summary.tsv", sep="\t", index=False)
    pd.DataFrame({
        "field": [
            "dataset", "author_annotated_release_cells", "training_control_cells", "heldout_perturbation_cells",
            "timepoint_networks", "network_pooling", "score_pooling", "cell_cap",
            "gene_cap", "region_cap", "native_output", "cross_timepoint_score_comparison",
            "causal_evaluation_scope", "scenicplus_cores",
            "offline_genome_manifest_sha256", "method_native_region_set_selection",
        ],
        "value": [
            "GSE274113", 137604, 11447, 126157, 4, "none", "none", "none", "none",
            "none", "independent pycisTopic and SCENIC+ outputs per timepoint",
            "prohibited", "replicate-paired heldout effects evaluated within TF-timepoint strata",
            int(settings['scenicplus_cores']),
            sha256(Path(settings['offline_genome_manifest'])),
            "Otsu and top-3000 regions per selected topic; both collections retained",
        ],
    }).to_csv(root / "run_validation.tsv", sep="\t", index=False)


    rows = []
    for path in sorted(p for p in root.rglob("*") if p.is_file() and "tmp" not in p.parts):
        if path == root / "artifact_manifest.tsv":
            continue
        rows.append({"file": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(rows).to_csv(root / "artifact_manifest.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


def run_scenicplus(object_root=OBJECT_ROOT, output_parent=OUTPUT_PARENT,
                   log_root=LOG_ROOT, offline_producer=OFFLINE_PRODUCER,
                   cores=CORES, monitor=True):
    if cores < 1:
        raise ValueError("cores must be positive")
    object_root, output_parent, log_root = map(Path, (object_root, output_parent, log_root))
    source_root = EXECUTION_ROOT / "sources/scenicplus"
    resources = REPO / "data/gse274113/external_resources/SCENICplus_hg38_screen_v10_clust"
    offline = REPO / "data/gse274113/external_resources/SCENICplus_GRCh38p14_Ensembl116"
    destination = output_parent / "SCENICplus"
    work = output_parent / ".SCENICplus.work"
    entry_script = Path(__file__).absolute()
    settings = {
        "ranking_db": str(resources / "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"),
        "scores_db": str(resources / "hg38_screen_v10_clust.regions_vs_motifs.scores.feather"),
        "motif_annotations": str(resources / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"),
        "resource_manifest": str(resources / "artifact_manifest.tsv"),
        "offline_genome_annotation": str(offline / "genome_annotation.tsv"),
        "offline_chromsizes": str(offline / "chromsizes.tsv"),
        "offline_genome_manifest": str(offline / "artifact_manifest.tsv"),
        "offline_resource_root": str(offline), "offline_resource_producer": str(offline_producer),
        "snakefile": str(source_root / "src/scenicplus/snakemake/Snakefile"),
        "cfg_default": str(source_root / "src/scenicplus/snakemake/config.yaml"),
        "default_config": str(source_root / "src/scenicplus/snakemake/config.yaml"),
        "source_commit_file": str(REPO / "logs/gse274113/comparator_install/scenicplus_source_commit.txt"),
        "producer_script": str(entry_script), "scenicplus_cores": cores,
    }
    required = [Path(settings[name]) for name in (
        "ranking_db", "scores_db", "motif_annotations", "resource_manifest",
        "offline_genome_annotation", "offline_chromsizes", "offline_genome_manifest",
        "offline_resource_producer", "snakefile", "default_config", "source_commit_file")]
    required += [offline / "resource_validation.tsv", offline / "source_manifest.tsv",
                 object_root / "input/input_validation.tsv", object_root / "cistopic_models/run_validation.tsv"]
    for path in required:
        if not path.is_file() or not path.stat().st_size:
            raise RuntimeError(f"Missing formal SCENIC+ input: {path}")
    validate_offline_genome_resources(settings)
    recorded = Path(settings["source_commit_file"]).read_text().strip()
    current = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source_root, text=True).strip()
    dirty = subprocess.check_output(["git", "status", "--porcelain"], cwd=source_root, text=True).strip()
    if recorded != current or dirty:
        raise RuntimeError("SCENIC+ source tree is not clean at the recorded commit")
    if destination.exists():
        raise RuntimeError(f"Formal SCENIC+ result already exists: {destination}")
    work.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    python = python_executable("scenicplus")
    env = task_environment(python=python)
    with tee_log(log_root / f"run_scenicplus_gse274113_timepoint_full_{stamp}.log"):
        print(f"started_utc\t{utc_now()}")
        print("network_unit\tfour independently fitted control-timepoint eRegulon networks\n"
              "network_pooling\tnone\ncell_cap\tnone\ngene_cap\tnone\nregion_cap\tnone")
        for index, (slug, label, cells) in enumerate(TIMEPOINTS):
            time_work = work / slug
            settings.update({
                "time_work": str(time_work), "timepoint": label, "expected_cells": cells,
                "input_cistopic": str(object_root / f"cistopic_models/{slug}/GSE274113_NT_{slug}_selected_cistopic.pkl"),
                "input_gex": str(object_root / f"input/{slug}/GSE274113_NT_{slug}_GEX.h5ad"),
                "region_sets": str(object_root / f"cistopic_models/{slug}/region_sets"),
                "cfg_path": str(time_work / "config.yaml"), "work_tmp": str(time_work / "tmp"),
                "seed": 666 + index,
            })
            time_work.mkdir(parents=True, exist_ok=True)
            for name in ("input_cistopic", "input_gex", "region_sets"):
                if not Path(settings[name]).exists():
                    raise RuntimeError(f"Missing {label} SCENIC+ input: {settings[name]}")
            markers = [time_work / "run_validation.tsv", time_work / "artifact_manifest.tsv"]
            if all(path.is_file() and path.stat().st_size for path in markers):
                validate_timepoint_artifacts(settings)
                shutil.rmtree(time_work / "tmp", ignore_errors=True)
                print(f"existing_validated\t{label}")
                continue
            if any(path.exists() for path in markers):
                raise RuntimeError(f"Incomplete SCENIC+ resume markers for {label}")
            (time_work / "tmp").mkdir(exist_ok=True)
            for name, source in (("genome_annotation.tsv", offline / "genome_annotation.tsv"),
                                 ("chromsizes.tsv", offline / "chromsizes.tsv")):
                target = time_work / name
                if target.exists() and target.read_bytes() != source.read_bytes():
                    raise RuntimeError(f"Existing {name} differs from offline resource for {label}")
                if not target.exists():
                    shutil.copyfile(source, target)
                    target.chmod(0o444)
            write_config(settings)
            print(f"timepoint_started\t{label}\t{utc_now()}")
            started = int(time.time())
            run_task([Path(python).parent / "snakemake", "--snakefile", settings["snakefile"],
                      "--configfile", settings["cfg_path"], "--cores", str(cores),
                      "--rerun-incomplete", "--keep-going", "--latency-wait", "120", "--printshellcmds"],
                     env=env, cwd=time_work, stage=f"gse274113_scenicplus_eregulon_{slug}",
                     workers=cores, producer=entry_script, monitor=monitor, timed=True)
            settings["fit_seconds"] = int(time.time()) - started
            write_timepoint_records(settings)
            validate_timepoint_artifacts(settings)
            shutil.rmtree(time_work / "tmp")
            print(f"timepoint_completed\t{label}\t{utc_now()}")
        settings["work_dir"] = str(work)
        write_suite_records(settings)
        work.rename(destination)
        print(f"destination\t{destination}\ncompleted_utc\t{utc_now()}")


def main():
    run_scenicplus(**TASK_SETTINGS["network"])


def run_timing():
    run_scenicplus(**TASK_SETTINGS["timing"])


if __name__ == "__main__":
    use_python("scenicplus")
    main()
