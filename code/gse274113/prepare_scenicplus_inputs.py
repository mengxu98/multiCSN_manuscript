#!/usr/bin/env python3
from __future__ import annotations

import os
import pickle
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "functions"))
from runtime import use_python
from hdf5_io import read_sparse_group
OBJECT_ROOT = Path(__file__).absolute().parents[2] / "data/gse274113/objects/GSE274113_SCENICplus_full"
MALLET_MEMORY = "120G"
if __name__ == "__main__":
    use_python("scenicplus", mallet_memory=MALLET_MEMORY)

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from pycisTopic.cistopic_class import create_cistopic_object


TIMEPOINTS = {
    "day 7": 3304,
    "day 9": 4167,
    "day 11": 1903,
    "day 14": 2073,
}


def canonicalize_region_names(region_names: list[str]) -> list[str]:
    canonical: list[str] = []
    for region in region_names:
        if ":" in region:
            chromosome, interval = region.split(":", 1)
            try:
                start_text, end_text = interval.split("-", 1)
            except ValueError as error:
                raise RuntimeError(f"Invalid region identifier: {region}") from error
        else:
            try:
                chromosome, start_text, end_text = region.rsplit("-", 2)
            except ValueError as error:
                raise RuntimeError(f"Invalid region identifier: {region}") from error
        try:
            start = int(start_text)
            end = int(end_text)
        except ValueError as error:
            raise RuntimeError(f"Non-integer region coordinates: {region}") from error
        if not chromosome or start < 0 or end <= start:
            raise RuntimeError(f"Invalid region coordinates: {region}")
        canonical.append(f"{chromosome}:{start}-{end}")
    if len(set(canonical)) != len(canonical):
        raise RuntimeError("Region canonicalization produced duplicated identifiers.")
    return canonical


def main() -> None:
    script = Path(__file__).resolve()
    repo_root = script.parents[2]
    object_root = repo_root / "data/gse274113/objects/GSE274113_control_full"
    counts_path = object_root / "GSE274113_all_NT_controls_counts.h5"
    metadata_path = object_root / "GSE274113_all_NT_controls_metadata.tsv.gz"
    validation_path = object_root / "control_object_validation.tsv"
    blacklist_path = (
        repo_root
        / "data/gse274113/external_resources/SCENICplus_auxiliary/hg38-blacklist.v2.bed"
    )
    parent = OBJECT_ROOT
    destination = parent / "input"
    required = [counts_path, metadata_path, validation_path, blacklist_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("Missing formal SCENIC+ inputs:\n" + "\n".join(missing))
    if destination.exists():
        raise SystemExit(f"Formal SCENIC+ input already exists: {destination}")

    parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".input.staging_", dir=parent))
    try:
        with h5py.File(counts_path, "r") as handle:
            rna, genes, rna_cells = read_sparse_group(handle, "rna")
            atac, regions, atac_cells = read_sparse_group(handle, "atac")

        source_validation = pd.read_csv(validation_path, sep="\t", dtype=str)
        if set(source_validation.columns) != {"field", "value"}:
            raise RuntimeError("Control-object validation schema is invalid.")
        source_validation = source_validation.set_index("field")["value"]
        expected_rna_features = int(source_validation["rna_release_features"])
        expected_atac_features = int(source_validation["atac_consensus_features"])
        expected_cells = int(source_validation["network_training_control_cells"])
        if rna.shape != (expected_rna_features, expected_cells) or atac.shape != (
            expected_atac_features, expected_cells
        ):
            raise RuntimeError(f"Unexpected full matrices: RNA={rna.shape}, ATAC={atac.shape}")
        if rna_cells != atac_cells or len(set(rna_cells)) != expected_cells:
            raise RuntimeError("RNA and ATAC cell identities are not an exact one-to-one match.")
        if len(set(genes)) != len(genes) or len(set(regions)) != len(regions):
            raise RuntimeError("SCENIC+ input feature identifiers are not unique.")
        if np.any(rna.data < 0) or np.any(atac.data < 0):
            raise RuntimeError("Count matrices contain negative entries.")
        regions = canonicalize_region_names(regions)

        metadata = pd.read_csv(metadata_path, sep="\t", low_memory=False)
        if "cell_id" not in metadata or metadata["cell_id"].duplicated().any():
            raise RuntimeError("Exported control metadata lacks unique cell identifiers.")
        metadata = metadata.set_index("cell_id", drop=False).loc[rna_cells].copy()
        if (
            len(metadata) != 11447
            or set(metadata["Timepoint"]) != {"day 7", "day 9", "day 11", "day 14"}
            or metadata["replicate"].nunique() != 14
            or set(metadata["target"]) != {"NT"}
        ):
            raise RuntimeError("Exported metadata violates the complete control-cell contract.")

        cell_array = np.asarray(rna_cells, dtype=object)
        summary_rows: list[dict[str, object]] = []
        for timepoint, expected_cells in TIMEPOINTS.items():
            slug = timepoint.replace(" ", "_")
            time_dir = stage / slug
            time_dir.mkdir()
            mask = metadata["Timepoint"].to_numpy() == timepoint
            time_metadata = metadata.loc[mask].copy()
            time_cells = cell_array[mask].tolist()
            if len(time_cells) != expected_cells or time_metadata["replicate"].nunique() < 1:
                raise RuntimeError(f"Unexpected full control input for {timepoint}.")

            time_rna = rna[:, mask].transpose().tocsr()
            time_atac = atac[:, mask].tocsr()
            gex = ad.AnnData(
                X=time_rna,
                obs=time_metadata,
                var=pd.DataFrame(index=pd.Index(genes, name="gene")),
            )
            gex.raw = gex.copy()
            sc.pp.normalize_total(gex, target_sum=1e4)
            sc.pp.log1p(gex)
            gex_path = time_dir / f"GSE274113_NT_{slug}_GEX.h5ad"
            gex.write_h5ad(gex_path, compression="gzip")

            cistopic = create_cistopic_object(
                fragment_matrix=time_atac,
                cell_names=time_cells,
                region_names=regions,
                path_to_blacklist=str(blacklist_path),
                min_frag=1,
                min_cell=1,
                is_acc=1,
                project=f"GSE274113_NT_{slug}",
                tag_cells=False,
            )
            cistopic.add_cell_data(time_metadata, split_pattern="___")
            if len(cistopic.cell_names) != expected_cells or set(cistopic.cell_names) != set(time_cells):
                raise RuntimeError(f"pycisTopic preprocessing lost cells at {timepoint}.")
            if len(cistopic.region_names) <= 0 or len(cistopic.region_names) > len(regions):
                raise RuntimeError(f"pycisTopic preprocessing produced invalid regions at {timepoint}.")
            cistopic_path = time_dir / f"GSE274113_NT_{slug}_cistopic.pkl"
            with cistopic_path.open("wb") as handle:
                pickle.dump(cistopic, handle, protocol=pickle.HIGHEST_PROTOCOL)
            summary_rows.append(
                {
                    "timepoint": timepoint,
                    "control_cells": expected_cells,
                    "replicates": time_metadata["replicate"].nunique(),
                    "input_rna_features": len(genes),
                    "input_atac_consensus_features": len(regions),
                    "pycistopic_regions_after_native_filtering": len(cistopic.region_names),
                    "pycistopic_binary_nonzero": cistopic.binary_matrix.nnz,
                    "network_pooling": "none",
                    "cell_cap": "none",
                    "gene_cap": "none",
                    "region_cap": "none",
                    "pycistopic_input_filter": (
                        "hg38 blacklist exclusion plus pycisTopic min_frag=1 and min_cell=1"
                    ),
                    "region_identifier_normalization": (
                        "separator-only canonicalization to chrom:start-end; coordinates and order unchanged"
                    ),
                }
            )
            del gex, cistopic, time_rna, time_atac

        summary = pd.DataFrame(summary_rows)
        if list(summary["timepoint"]) != list(TIMEPOINTS) or int(summary["control_cells"].sum()) != 11447:
            raise RuntimeError("Independent SCENIC+ input suite failed full-cell accounting.")
        summary.to_csv(stage / "timepoint_input_summary.tsv", sep="\t", index=False)

        os.replace(stage, destination)
        print(f"destination\t{destination}")
    except Exception:
        if stage.exists():
            shutil.rmtree(stage)
        raise


if __name__ == "__main__":
    main()
