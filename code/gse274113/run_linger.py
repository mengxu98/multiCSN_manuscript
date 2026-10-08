#!/usr/bin/env python3
"""Run the LINGER 1.110 workflow on one GSE274113 control time point."""

from __future__ import annotations

import csv
import gzip
import importlib.metadata
import json
import os
import random
import sys
import time
import types
from pathlib import Path

# Invoked per timepoint (LingerGRN 1.110, seed 42).
ARGUMENTS = {
    "counts": Path(
        "data/gse274113/objects/GSE274113_control_full/"
        "GSE274113_all_NT_controls_counts.h5"
    ),
    "metadata": Path(
        "data/gse274113/objects/GSE274113_control_full/"
        "GSE274113_all_NT_controls_metadata.tsv.gz"
    ),
    "grn_dir": Path("data/resources/LINGER/data_bulk"),
    "output_root": Path("runs/gse274113/LINGER"),
    "timepoint": "day 7",  # day 7 / day 9 / day 11 / day 14
    "phase": "all",
    "seed": 42,
}


sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "functions"))
from runtime import use_python, sha256
from hdf5_io import read_sparse_group
if __name__ == "__main__":
    use_python("linger", hash_seed=ARGUMENTS["seed"])

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


EXPECTED_CELLS = {"day_7": 3304, "day_9": 4167, "day_11": 1903, "day_14": 2073}
PHASES = ("prepare", "preprocess", "train", "networks", "export")


def canonical_region(region: str) -> str:
    if ":" in region:
        chromosome, interval = region.split(":", 1)
        start, end = interval.split("-", 1)
    else:
        chromosome, start, end = region.rsplit("-", 2)
    start_value, end_value = int(start), int(end)
    if not chromosome or start_value < 0 or end_value <= start_value:
        raise ValueError(f"Invalid region: {region}")
    return f"{chromosome}:{start_value}-{end_value}"


def marker(root: Path, phase: str) -> Path:
    return root / f".{phase}.done.json"


def mark(root: Path, phase: str, started: float, details: dict[str, object]) -> None:
    marker(root, phase).write_text(
        json.dumps(
            {
                "phase": phase,
                "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(started)),
                "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "elapsed_seconds": round(time.time() - started, 3),
                **details,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def prepare(
    root: Path, counts_path: Path, metadata_path: Path, timepoint: str, seed: int
) -> None:
    if marker(root, "prepare").is_file():
        print(f"prepare already complete: {root}", flush=True)
        return
    started = time.time()
    data_dir = root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    with h5py.File(counts_path, "r") as handle:
        rna, genes, rna_cells = read_sparse_group(handle, "rna")
        atac, regions, atac_cells = read_sparse_group(handle, "atac")
    if rna.shape != (36601, 11447) or atac.shape != (226206, 11447):
        raise RuntimeError(f"Unexpected full matrices: RNA={rna.shape}, ATAC={atac.shape}")
    if rna_cells != atac_cells or len(set(rna_cells)) != 11447:
        raise RuntimeError("RNA and ATAC cell identities are not an exact one-to-one match")
    regions = [canonical_region(region) for region in regions]
    if len(set(genes)) != len(genes) or len(set(regions)) != len(regions):
        raise RuntimeError("Feature identifiers are not unique")

    metadata = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    if metadata["cell_id"].duplicated().any() or set(metadata["cell_id"]) != set(rna_cells):
        raise RuntimeError("Metadata does not match the control matrices")
    metadata = metadata.set_index("cell_id", drop=False).loc[rna_cells]
    expected_design = {key.replace("_", " "): value for key, value in EXPECTED_CELLS.items()}
    observed_design = metadata.groupby("Timepoint", observed=True).size().to_dict()
    if observed_design != expected_design or set(metadata["target"]) != {"NT"}:
        raise RuntimeError(f"Control-cell design changed: {observed_design}")

    label = timepoint.replace("_", " ")
    mask = metadata["Timepoint"].to_numpy() == label
    current_metadata = metadata.loc[mask].copy()
    current_metadata["barcode"] = current_metadata["cell_id"]
    current_metadata["sample"] = current_metadata["replicate"].astype(str)
    current_metadata["label"] = current_metadata["annotation_simplified"].fillna("Unannotated")
    if len(current_metadata) != EXPECTED_CELLS[timepoint]:
        raise RuntimeError(f"Unexpected cell count for {timepoint}")

    rna_data = ad.AnnData(
        X=rna[:, mask].T.tocsr(),
        obs=current_metadata,
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
    )
    atac_data = ad.AnnData(
        X=atac[:, mask].T.tocsr(),
        obs=current_metadata.copy(),
        var=pd.DataFrame(index=pd.Index(regions, name="region")),
    )
    rna_data.var["gene_ids"] = rna_data.var_names
    atac_data.var["gene_ids"] = atac_data.var_names

    sc.pp.filter_cells(rna_data, min_genes=200)
    sc.pp.filter_genes(rna_data, min_cells=3)
    sc.pp.filter_cells(atac_data, min_genes=200)
    sc.pp.filter_genes(atac_data, min_cells=3)
    common = [cell for cell in current_metadata.index if cell in rna_data.obs_names and cell in atac_data.obs_names]
    rna_data = rna_data[common].copy()
    atac_data = atac_data[common].copy()
    if rna_data.n_obs != atac_data.n_obs or not np.array_equal(rna_data.obs_names, atac_data.obs_names):
        raise RuntimeError("Modality filtering broke cell alignment")

    rna_data.write_h5ad(data_dir / "adata_RNA.h5ad", compression="gzip")
    atac_data.write_h5ad(data_dir / "adata_ATAC.h5ad", compression="gzip")

    from LingerGRN.pseudo_bulk import pseudo_bulk

    random.seed(seed)
    np.random.seed(seed)
    samples = sorted(atac_data.obs["sample"].unique())
    single_pseudobulk = len(samples) ** 2 > 100
    tg_parts: list[pd.DataFrame] = []
    re_parts: list[pd.DataFrame] = []
    for sample in samples:
        selected = rna_data.obs["sample"] == sample
        tg_part, re_part = pseudo_bulk(
            rna_data[selected].copy(), atac_data[selected].copy(), single_pseudobulk
        )
        tg_parts.append(tg_part)
        re_parts.append(re_part)
    tg_pseudobulk = pd.concat(tg_parts, axis=1).fillna(0)
    re_pseudobulk = pd.concat(re_parts, axis=1).fillna(0)
    re_pseudobulk[re_pseudobulk > 100] = 100
    if tg_pseudobulk.columns.duplicated().any() or re_pseudobulk.columns.duplicated().any():
        raise RuntimeError("Pseudobulk identifiers are duplicated")
    if list(tg_pseudobulk.columns) != list(re_pseudobulk.columns):
        raise RuntimeError("RNA and ATAC pseudobulk columns differ")
    if re_pseudobulk.index.duplicated().any():
        raise RuntimeError("ATAC pseudobulk regions are duplicated")
    pd.Series(re_pseudobulk.index).to_csv(data_dir / "Peaks.txt", sep="\t", header=False, index=False)
    tg_pseudobulk.to_csv(data_dir / "TG_pseudobulk.tsv", sep="\t")
    re_pseudobulk.to_csv(data_dir / "RE_pseudobulk.tsv", sep="\t")
    pd.DataFrame(
        {
            "field": [
                "dataset",
                "timepoint",
                "source_cells",
                "retained_cells",
                "replicates",
                "cell_type_labels",
                "rna_features",
                "atac_features",
                "pseudobulks",
                "seed",
            ],
            "value": [
                "GSE274113",
                label,
                EXPECTED_CELLS[timepoint],
                rna_data.n_obs,
                len(samples),
                rna_data.obs["label"].nunique(),
                rna_data.n_vars,
                atac_data.n_vars,
                tg_pseudobulk.shape[1],
                seed,
            ],
        }
    ).to_csv(root / "input_summary.tsv", sep="\t", index=False)
    mark(
        root,
        "prepare",
        started,
        {
            "counts_sha256": sha256(counts_path),
            "metadata_sha256": sha256(metadata_path),
            "retained_cells": rna_data.n_obs,
            "pseudobulks": tg_pseudobulk.shape[1],
        },
    )


def preprocess(root: Path, grn_dir: Path) -> None:
    if marker(root, "preprocess").is_file():
        print(f"preprocess already complete: {root}", flush=True)
        return
    if not marker(root, "prepare").is_file():
        raise RuntimeError("prepare phase is incomplete")
    started = time.time()
    data_dir, output_dir = root / "data", root / "output"
    output_dir.mkdir(exist_ok=True)
    tg = pd.read_csv(data_dir / "TG_pseudobulk.tsv", sep="\t", index_col=0)
    re = pd.read_csv(data_dir / "RE_pseudobulk.tsv", sep="\t", index_col=0)
    if tg.empty or re.empty or list(tg.columns) != list(re.columns):
        raise RuntimeError(f"Malformed or unaligned pseudobulk matrices: RNA={tg.shape}, ATAC={re.shape}")
    peaks = pd.read_csv(data_dir / "Peaks.txt", sep="\t", header=None)[0].astype(str)
    if list(peaks) != list(re.index.astype(str)):
        raise RuntimeError("Peaks.txt does not exactly match the ATAC pseudobulk index")
    from LingerGRN.preprocess import preprocess as linger_preprocess

    linger_preprocess(tg, re, f"{grn_dir}/", "hg38", "LINGER", f"{output_dir}/")
    required = ["Exp.txt", "Symbol.txt", "TFexp.txt", "TFName.txt", "Openness.txt", "TF_binding.txt", "index.txt"]
    missing = [name for name in required if not (output_dir / name).is_file()]
    if missing:
        raise RuntimeError(f"LINGER preprocessing outputs are missing: {missing}")
    mark(root, "preprocess", started, {"required_outputs": required})


def validate_training_outputs(root: Path) -> pd.DataFrame:
    import torch

    output_dir = root / "output"
    data_merge = pd.read_csv(output_dir / "data_merge.txt", sep="\t", index_col=0)
    if data_merge["Symbol"].duplicated().any():
        raise RuntimeError("LINGER data_merge contains duplicated gene symbols")
    rows = []
    for chromosome in [f"chr{index}" for index in range(1, 23)] + ["chrX"]:
        expected = set(map(int, data_merge.index[data_merge["chr"] == chromosome]))
        net = set(map(int, torch.load(output_dir / f"net_{chromosome}.pt", map_location="cpu", weights_only=False)))
        shap = set(map(int, torch.load(output_dir / f"shap_{chromosome}.pt", map_location="cpu", weights_only=False)))
        result = pd.read_csv(output_dir / f"result_{chromosome}.txt", sep="\t", index_col=0)
        failed = set(result.index[pd.to_numeric(result.iloc[:, 0], errors="coerce") == -100])
        skipped = expected - net
        skipped_symbols = {str(data_merge.loc[index, "Symbol"]) for index in skipped}
        if net != shap or net - expected or skipped_symbols != failed:
            raise RuntimeError(
                f"Inconsistent LINGER training keys for {chromosome}: "
                f"expected={len(expected)}, net={len(net)}, shap={len(shap)}, "
                f"skipped={sorted(skipped_symbols)}, failed={sorted(failed)}"
            )
        rows.append(
            {
                "chromosome": chromosome,
                "expected_genes": len(expected),
                "trained_genes": len(net),
                "method_skipped_genes": len(skipped),
                "method_skipped_symbols": ",".join(sorted(skipped_symbols)),
            }
        )
    coverage = pd.DataFrame(rows)
    return coverage


def train(root: Path, grn_dir: Path) -> None:
    if marker(root, "train").is_file():
        validate_training_outputs(root)
        print(f"train already complete: {root}", flush=True)
        return
    if not marker(root, "preprocess").is_file():
        raise RuntimeError("preprocess phase is incomplete")
    started = time.time()
    from LingerGRN import LINGER_tr

    output_dir = root / "output"
    LINGER_tr.training(f"{grn_dir}/", "LINGER", f"{output_dir}/", "ReLU", "Human")
    chromosomes = [f"chr{index}" for index in range(1, 23)] + ["chrX"]
    required = [output_dir / f"net_{chromosome}.pt" for chromosome in chromosomes]
    required += [output_dir / f"shap_{chromosome}.pt" for chromosome in chromosomes]
    missing = [path.name for path in required if not path.is_file()]
    if missing:
        raise RuntimeError(f"LINGER training outputs are missing: {missing}")
    coverage = validate_training_outputs(root)
    mark(
        root,
        "train",
        started,
        {
            "chromosomes": len(chromosomes),
            "trained_genes": int(coverage["trained_genes"].sum()),
            "method_skipped_genes": int(coverage["method_skipped_genes"].sum()),
        },
    )


def networks(root: Path, grn_dir: Path) -> None:
    if marker(root, "networks").is_file():
        print(f"networks already complete: {root}", flush=True)
        return
    if not marker(root, "train").is_file():
        raise RuntimeError("train phase is incomplete")
    started = time.time()
    data_dir, output_dir = root / "data", root / "output"
    rna_data = ad.read_h5ad(data_dir / "adata_RNA.h5ad")
    atac_data = ad.read_h5ad(data_dir / "adata_ATAC.h5ad")
    from LingerGRN import LL_net

    LL_net.TF_RE_binding(f"{grn_dir}/", rna_data, atac_data, "hg38", "LINGER", f"{output_dir}/")
    LL_net.cis_reg(f"{grn_dir}/", rna_data, atac_data, "hg38", "LINGER", f"{output_dir}/")
    LL_net.trans_reg(f"{grn_dir}/", "LINGER", f"{output_dir}/", "hg38")
    required = [
        "cell_population_TF_RE_binding.txt",
        "cell_population_cis_regulatory.txt",
        "cell_population_trans_regulatory.txt",
    ]
    missing = [name for name in required if not (output_dir / name).is_file()]
    if missing:
        raise RuntimeError(f"LINGER network outputs are missing: {missing}")
    mark(root, "networks", started, {"required_outputs": required})


def write_dense_edges(
    source: Path, destination: Path, row_name: str, column_name: str
) -> int:
    rows_written = 0
    with gzip.open(destination, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([row_name, column_name, "score"])
        for chunk in pd.read_csv(source, sep="\t", index_col=0, chunksize=500):
            stacked = chunk.stack(dropna=True)
            for (row, column), score in stacked.items():
                value = float(score)
                if np.isfinite(value) and value != 0:
                    writer.writerow([row, column, f"{value:.17g}"])
                    rows_written += 1
    return rows_written


def export(root: Path) -> None:
    if marker(root, "export").is_file():
        print(f"export already complete: {root}", flush=True)
        return
    if not marker(root, "networks").is_file():
        raise RuntimeError("networks phase is incomplete")
    started = time.time()
    output_dir, network_dir = root / "output", root / "networks"
    network_dir.mkdir(exist_ok=True)
    tf_region_rows = write_dense_edges(
        output_dir / "cell_population_TF_RE_binding.txt",
        network_dir / "tf_region.tsv.gz",
        "region",
        "regulator",
    )
    tf_gene_rows = write_dense_edges(
        output_dir / "cell_population_trans_regulatory.txt",
        network_dir / "tf_gene.tsv.gz",
        "target",
        "regulator",
    )
    region_gene_rows = 0
    with gzip.open(network_dir / "region_gene.tsv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["region", "target", "score"])
        source = output_dir / "cell_population_cis_regulatory.txt"
        for chunk in pd.read_csv(source, sep="\t", header=None, names=["region", "target", "score"], chunksize=500_000):
            chunk["score"] = pd.to_numeric(chunk["score"], errors="coerce")
            chunk = chunk[np.isfinite(chunk["score"]) & (chunk["score"] != 0)]
            for row in chunk.itertuples(index=False):
                writer.writerow([row.region, row.target, f"{float(row.score):.17g}"])
                region_gene_rows += 1
    summary = pd.DataFrame(
        {
            "endpoint": ["tf_gene", "tf_region", "region_gene"],
            "relations": [tf_gene_rows, tf_region_rows, region_gene_rows],
        }
    )
    summary.to_csv(network_dir / "network_summary.tsv", sep="\t", index=False)
    mark(
        root,
        "export",
        started,
        {
            "tf_gene_relations": tf_gene_rows,
            "tf_region_relations": tf_region_rows,
            "region_gene_relations": region_gene_rows,
        },
    )


# ---- Parameters ----


def main(settings=None) -> None:
    settings = dict(ARGUMENTS if settings is None else settings)
    args = types.SimpleNamespace(**settings)
    # Resolve every path before the workflow changes into its output directory;
    # otherwise the relative in-script defaults would no longer resolve.
    for name in ("counts", "metadata", "grn_dir", "output_root"):
        setattr(args, name, Path(getattr(args, name)).resolve())
    if args.timepoint not in EXPECTED_CELLS:
        raise SystemExit(
            f"timepoint must be one of {sorted(EXPECTED_CELLS)}, got {args.timepoint}"
        )
    installed = importlib.metadata.version("LingerGRN")
    if installed != "1.110":
        raise SystemExit(f"Expected LingerGRN 1.110, found {installed}")
    for path in (args.counts, args.metadata):
        if not path.is_file():
            raise SystemExit(f"Missing input: {path}")
    if args.phase != "prepare" and not args.grn_dir.is_dir():
        raise SystemExit(f"Missing LINGER data_bulk directory: {args.grn_dir}")
    root = (args.output_root / args.timepoint.replace(" ", "_")).resolve()
    root.mkdir(parents=True, exist_ok=True)
    os.chdir(root)
    functions = {
        "prepare": lambda: prepare(root, args.counts, args.metadata, args.timepoint, args.seed),
        "preprocess": lambda: preprocess(root, args.grn_dir),
        "train": lambda: train(root, args.grn_dir),
        "networks": lambda: networks(root, args.grn_dir),
        "export": lambda: export(root),
    }
    selected = PHASES if args.phase == "all" else (args.phase,)
    for phase in selected:
        print(f"phase_started\t{phase}\t{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}", flush=True)
        functions[phase]()
        print(f"phase_finished\t{phase}\t{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}", flush=True)


def run_timepoint(timepoint, output_root):
    main(dict(ARGUMENTS, timepoint=timepoint, output_root=Path(output_root)))


if __name__ == "__main__":
    main()
