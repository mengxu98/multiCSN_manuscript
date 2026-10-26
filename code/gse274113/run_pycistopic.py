#!/usr/bin/env python3
from __future__ import annotations

import os
import pickle
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "functions"))
from runtime import use_python
OBJECT_ROOT = Path(__file__).absolute().parents[2] / "data/gse274113/objects/GSE274113_SCENICplus_full"
MALLET_MEMORY = "120G"
if __name__ == "__main__":
    use_python("scenicplus", mallet_memory=MALLET_MEMORY)

import pandas as pd
from pycisTopic.lda_models import evaluate_models, run_cgs_models_mallet
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.utils import region_names_to_coordinates


TOPICS = [2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
TIMEPOINTS = ("day 7", "day 9", "day 11", "day 14")


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def write_beds(collection: dict, directory: Path) -> list[dict]:
    directory.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for name, table in collection.items():
        regions = list(table.index)
        coordinates = region_names_to_coordinates(regions).sort_values(
            ["Chromosome", "Start", "End"]
        )
        path = directory / f"{safe_name(name)}.bed"
        coordinates.to_csv(path, sep="\t", header=False, index=False)
        rows.append(
            {
                "collection": directory.name,
                "name": str(name),
                "regions": len(regions),
                "file": str(path),
            }
        )
    return rows


def model_metrics(models: list) -> pd.DataFrame:
    rows = []
    for model in sorted(models, key=lambda item: item.n_topic):
        row = {"n_topics": model.n_topic}
        for column in model.metrics.columns:
            row[column] = model.metrics.loc["Metric", column]
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    script = Path(__file__).resolve()
    repo_root = script.parents[2]
    object_root = OBJECT_ROOT
    input_root = object_root / "input"
    metadata_path = repo_root / "data/gse274113/objects/GSE274113_control_full/GSE274113_all_NT_controls_metadata.tsv.gz"
    auxiliary = repo_root / "data/gse274113/external_resources/SCENICplus_auxiliary"
    mallet_path = auxiliary / "Mallet-202108/bin/mallet"
    parent = object_root
    destination = parent / "cistopic_models"
    work = parent / ".cistopic_models.work"
    input_paths = {
        timepoint: input_root / timepoint.replace(" ", "_")
        / f"GSE274113_NT_{timepoint.replace(' ', '_')}_cistopic.pkl"
        for timepoint in TIMEPOINTS
    }
    required = [mallet_path, metadata_path, *input_paths.values()]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Missing formal pycisTopic inputs:\n" + "\n".join(missing))
    if destination.exists():
        raise SystemExit(f"Formal pycisTopic result already exists: {destination}")
    metadata = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    columns = ["cell_id", "Timepoint", "target"]
    if metadata[columns].isna().any().any() or metadata["cell_id"].duplicated().any():
        raise RuntimeError("Missing or duplicate control-cell IDs.")
    if set(metadata["Timepoint"]) != set(TIMEPOINTS) or set(metadata["target"]) != {"NT"}:
        raise RuntimeError("Control metadata violates the timepoint design.")

    work.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, object]] = []
    for index, timepoint in enumerate(TIMEPOINTS):
        timepoint_started = time.monotonic()
        slug = timepoint.replace(" ", "_")
        time_work = work / slug
        time_work.mkdir(exist_ok=True)
        model_dir = time_work / "models"
        mallet_tmp = time_work / "mallet_tmp"
        model_dir.mkdir(exist_ok=True)
        mallet_tmp.mkdir(exist_ok=True)
        input_path = input_paths[timepoint]
        with input_path.open("rb") as handle:
            cistopic = pickle.load(handle)
        expected_cells = set(metadata.loc[metadata["Timepoint"] == timepoint, "cell_id"])
        if len(cistopic.cell_names) != len(expected_cells) or set(cistopic.cell_names) != expected_cells or set(cistopic.cell_data["Timepoint"]) != {timepoint}:
            raise RuntimeError(f"Invalid pycisTopic input for {timepoint}.")

        models = []
        for n_topics in TOPICS:
            model_path = model_dir / f"Topic{n_topics}.pkl"
            if model_path.is_file():
                with model_path.open("rb") as handle:
                    model = pickle.load(handle)
            else:
                reuse_corpus = (mallet_tmp / "corpus.mallet").is_file()
                model = run_cgs_models_mallet(
                    cistopic,
                    n_topics=[n_topics],
                    n_cpu=1,
                    n_iter=500,
                    random_state=555 + index,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    top_topics_coh=5,
                    tmp_path=str(mallet_tmp),
                    save_path=str(model_dir),
                    reuse_corpus=reuse_corpus,
                    mallet_path=str(mallet_path),
                )[0]
            if model.n_topic != n_topics:
                raise RuntimeError(f"Topic model mismatch for {timepoint}, {n_topics} topics.")
            models.append(model)

        metrics = model_metrics(models)
        metrics.to_csv(time_work / "model_metrics.tsv", sep="\t", index=False)
        selected = evaluate_models(
            models,
            select_model=None,
            return_model=True,
            metrics=["Minmo_2011", "loglikelihood"],
            min_topics_coh=5,
            plot=False,
            save=str(time_work / "model_selection.pdf"),
        )
        if selected.n_topic not in TOPICS:
            raise RuntimeError(f"Model selection failed for {timepoint}.")
        cistopic.add_LDA_model(selected)
        selected_path = time_work / f"GSE274113_NT_{slug}_selected_cistopic.pkl"
        with selected_path.open("wb") as handle:
            pickle.dump(cistopic, handle, protocol=pickle.HIGHEST_PROTOCOL)

        topics_otsu = binarize_topics(
            cistopic, target="region", method="otsu", smooth_topics=True, plot=False
        )
        topics_top3k = binarize_topics(
            cistopic, target="region", method="ntop", smooth_topics=True,
            ntop=3000, plot=False,
        )
        region_root = time_work / "region_sets"
        region_rows = []
        region_rows.extend(write_beds(topics_otsu, region_root / "Topics_otsu"))
        region_rows.extend(write_beds(topics_top3k, region_root / "Topics_top_3k"))
        region_summary = pd.DataFrame(region_rows)
        if region_summary.empty or (region_summary["regions"] <= 0).any():
            raise RuntimeError(f"One or more pycisTopic region sets are empty for {timepoint}.")
        region_summary.to_csv(time_work / "region_set_summary.tsv", sep="\t", index=False)
        otsu_unique_regions = len(
            set().union(*(set(table.index) for table in topics_otsu.values()))
        )
        top3000_unique_regions = len(
            set().union(*(set(table.index) for table in topics_top3k.values()))
        )

        summary_rows.append(
            {
                "timepoint": timepoint,
                "training_cells": len(cistopic.cell_names),
                "input_regions_after_native_filtering": len(cistopic.region_names),
                "candidate_topic_models": len(TOPICS),
                "selected_topics": selected.n_topic,
                "topic_otsu_sets": len(topics_otsu),
                "topic_top3000_sets": len(topics_top3k),
                "topic_otsu_unique_regions": otsu_unique_regions,
                "topic_top3000_unique_regions": top3000_unique_regions,
                "method_native_region_set_selection": (
                    "Otsu and top-3000 regions per selected topic; both collections retained"
                ),
                "topic_model_and_region_set_seconds": time.monotonic() - timepoint_started,
                "network_pooling": "none",
                "cell_cap": "none",
                "region_cap": "none",
            }
        )
        del cistopic, models, selected, topics_otsu, topics_top3k

    summary = pd.DataFrame(summary_rows)
    if list(summary["timepoint"]) != list(TIMEPOINTS) or int(summary["training_cells"].sum()) != len(metadata):
        raise RuntimeError("pycisTopic cell totals do not match.")
    summary.to_csv(work / "timepoint_model_summary.tsv", sep="\t", index=False)

    os.replace(work, destination)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
