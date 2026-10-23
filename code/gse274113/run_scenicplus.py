#!/usr/bin/env python3
"""Run four independent SCENIC+ timepoint networks with validated offline inputs."""
from __future__ import annotations
from pathlib import Path
import shutil
import sys
import time
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "functions"))
from runtime import (REPO, EXECUTION_ROOT, python_executable, use_python,
                     run_task, task_environment, utc_now)

CORES = 12
TIMEPOINTS = [("day_7", "day 7", 3304), ("day_9", "day 9", 4167),
              ("day_11", "day 11", 1903), ("day_14", "day 14", 2073)]
OBJECT_ROOT = REPO / "data/gse274113/objects/GSE274113_SCENICplus_full"
OUTPUT_PARENT = REPO / "runs/gse274113/method_fits"
TIMING_ROOT = REPO / "runs/gse274113/timing_benchmark/SCENICplus"
TASK_SETTINGS = {
    "network": dict(object_root=OBJECT_ROOT, output_parent=OUTPUT_PARENT,
                    cores=CORES, monitor=True),
    "timing": dict(object_root=OBJECT_ROOT, output_parent=TIMING_ROOT / "output",
                   cores=CORES, monitor=False),
}

def validate_offline_genome_resources(settings):
    import pandas as pd

    annotation = pd.read_csv(settings["offline_genome_annotation"], sep="\t")
    columns = {"Chromosome", "Start", "End", "Strand", "Gene", "Transcription_Start_Site", "Transcript_type"}
    if not columns.issubset(annotation.columns) or annotation.empty:
        raise ValueError("Invalid genome annotation columns or empty annotation")
    chromsizes = pd.read_csv(settings["offline_chromsizes"], sep="\t")
    expected = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
    if set(chromsizes["Chromosome"]) != expected or chromsizes["Chromosome"].duplicated().any():
        raise ValueError("Chromosome sizes must cover the 25 primary chromosomes")
    if not ((chromsizes["Start"] == 0) & (chromsizes["End"] > 0)).all():
        raise ValueError("Invalid chromosome bounds")


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


def summarize_timepoint(settings):
    root = Path(settings["time_work"])
    counts = {}
    for name in ("tf_to_gene_adj.tsv", "region_to_gene_adj.tsv",
                 "eRegulon_direct.tsv", "eRegulons_extended.tsv"):
        with (root / name).open("rt", errors="replace") as handle:
            counts[name] = max(sum(1 for _ in handle) - 1, 0)
    if any(value <= 0 for value in counts.values()):
        raise RuntimeError(f"One or more SCENIC+ edge tables are empty: {counts}")
    return {
        "timepoint": settings["timepoint"],
        "training_cells": int(settings["expected_cells"]),
        "tf_to_gene_rows": counts["tf_to_gene_adj.tsv"],
        "region_to_gene_rows": counts["region_to_gene_adj.tsv"],
        "direct_eRegulon_rows": counts["eRegulon_direct.tsv"],
        "extended_eRegulon_rows": counts["eRegulons_extended.tsv"],
        "scenicplus_fit_seconds": float(settings["fit_seconds"]),
        "scenicplus_cores": int(settings["scenicplus_cores"]),
        "network_pooling": "none", "cell_cap": "none", "gene_cap": "none", "region_cap": "none",
        "method_native_region_set_selection": "Otsu and top-3000 regions per selected topic; both collections retained",
    }


def write_suite_summary(root, rows):
    import pandas as pd

    summary = pd.DataFrame(rows)
    if int(summary["training_cells"].sum()) != 11447 or len(summary) != 4:
        raise RuntimeError("SCENIC+ independent timepoint suite failed full-cell accounting.")
    summary.to_csv(root / "timepoint_run_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


def run_scenicplus(object_root=OBJECT_ROOT, output_parent=OUTPUT_PARENT,
                   cores=CORES, monitor=True):
    if cores < 1:
        raise ValueError("cores must be positive")
    object_root, output_parent = map(Path, (object_root, output_parent))
    source_root = EXECUTION_ROOT / "sources/scenicplus"
    resources = REPO / "data/gse274113/external_resources/SCENICplus_hg38_screen_v10_clust"
    offline = REPO / "data/gse274113/external_resources/SCENICplus_GRCh38p14_Ensembl116"
    destination = output_parent / "SCENICplus"
    work = output_parent / ".SCENICplus.work"
    settings = {
        "ranking_db": str(resources / "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"),
        "scores_db": str(resources / "hg38_screen_v10_clust.regions_vs_motifs.scores.feather"),
        "motif_annotations": str(resources / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"),
        "offline_genome_annotation": str(offline / "genome_annotation.tsv"),
        "offline_chromsizes": str(offline / "chromsizes.tsv"),
        "snakefile": str(source_root / "src/scenicplus/snakemake/Snakefile"),
        "cfg_default": str(source_root / "src/scenicplus/snakemake/config.yaml"),
        "scenicplus_cores": cores,
    }
    required = [Path(settings[name]) for name in (
        "ranking_db", "scores_db", "motif_annotations",
        "offline_genome_annotation", "offline_chromsizes", "snakefile", "cfg_default")]
    for path in required:
        if not path.is_file() or not path.stat().st_size:
            raise RuntimeError(f"Missing formal SCENIC+ input: {path}")
    validate_offline_genome_resources(settings)
    if destination.exists():
        raise RuntimeError(f"Formal SCENIC+ result already exists: {destination}")
    work.mkdir(parents=True, exist_ok=True)
    python = python_executable("scenicplus")
    env = task_environment(python=python)
    print(f"started_utc\t{utc_now()}")
    print("network_unit\tfour independently fitted control-timepoint eRegulon networks\n"
          "network_pooling\tnone\ncell_cap\tnone\ngene_cap\tnone\nregion_cap\tnone")
    rows = []
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
        (time_work / "tmp").mkdir(exist_ok=True)
        for name, source in (("genome_annotation.tsv", offline / "genome_annotation.tsv"),
                             ("chromsizes.tsv", offline / "chromsizes.tsv")):
            target = time_work / name
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
                 workers=cores, monitor=monitor, timed=True)
        settings["fit_seconds"] = int(time.time()) - started
        rows.append(summarize_timepoint(settings))
        shutil.rmtree(time_work / "tmp")
        print(f"timepoint_completed\t{label}\t{utc_now()}")
    write_suite_summary(work, rows)
    for slug, _, _ in TIMEPOINTS:
        workflow = work / slug / ".snakemake"
        if workflow.is_dir():
            shutil.rmtree(workflow)
    work.rename(destination)
    print(f"destination\t{destination}\ncompleted_utc\t{utc_now()}")


def main():
    run_scenicplus(**TASK_SETTINGS["network"])


def run_timing():
    run_scenicplus(**TASK_SETTINGS["timing"])


if __name__ == "__main__":
    use_python("scenicplus")
    main()
