from pathlib import Path
import tempfile
import gzip
import itertools
import warnings

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_count_chunks(path, *, header, index_col, columns, chunk_rows=512):
    blocks, features = [], []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        if header is not None:
            observed = next(handle).rstrip("\r\n").split("\t")
            if len(observed) == len(columns) + 1:
                observed = observed[1:]
            if observed != columns:
                raise ValueError(f"Count-matrix columns differ from the released paired cells: {path}")
        while True:
            lines = list(itertools.islice(handle, chunk_rows))
            if not lines:
                break
            if index_col is not None:
                pairs = [line.split("\t", 1) for line in lines]
                if any(len(pair) != 2 or not pair[0] for pair in pairs):
                    raise ValueError(f"Invalid released gene identifiers: {path}")
                features.extend(pair[0] for pair in pairs)
                numeric = "".join(pair[1] for pair in pairs)
            else:
                numeric = "".join(lines)
            with warnings.catch_warnings():
                warnings.simplefilter("error", DeprecationWarning)
                values = np.fromstring(numeric, sep="\t", dtype=np.int32)
            if values.size != len(lines) * len(columns) or (values < 0).any():
                raise ValueError(f"Invalid count-matrix dimensions or counts: {path}")
            blocks.append(sp.csr_matrix(values.reshape(len(lines), len(columns))))
    if not blocks:
        raise ValueError(f"Empty count matrix: {path}")
    return sp.vstack(blocks, format="csc"), features


def prepare_paired_counts(raw_root, output):
    raw_root, output = Path(raw_root), Path(output)
    if output.exists():
        raise FileExistsError(output)
    metadata = pd.read_csv(raw_root / "GSE162170_multiome_cell_metadata.txt.gz", sep="\t")
    cells = metadata["Cell.ID"].astype(str).tolist()
    if len(cells) != len(set(cells)):
        raise ValueError("Duplicate paired release cells")
    peaks = pd.read_csv(raw_root / "GSE162170_multiome_atac_consensus_peaks.txt.gz", sep="\t")
    peak_names = (peaks.seqnames.astype(str) + "-" + peaks.start.astype(str) + "-" + peaks.end.astype(str)).tolist()
    rna, genes = read_count_chunks(raw_root / "GSE162170_multiome_rna_counts.tsv.gz", header=0, index_col=0, columns=cells)
    atac, _ = read_count_chunks(raw_root / "GSE162170_multiome_atac_counts.tsv.gz", header=None, index_col=None, columns=cells)
    if atac.shape[0] != len(peak_names) or len(genes) != rna.shape[0] or len(set(peak_names)) != len(peak_names):
        raise ValueError("Released paired feature counts or peak coordinates differ")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".paired-counts-", dir=output.parent) as temporary:
        staged = Path(temporary) / "counts.h5"
        with h5py.File(staged, "w") as handle:
            for name, matrix, features in (("rna", rna, genes), ("atac", atac, peak_names)):
                group = handle.create_group(name)
                for key, values in (("data", matrix.data), ("indices", matrix.indices), ("indptr", matrix.indptr), ("shape", matrix.shape)):
                    group.create_dataset(key, data=values)
                group.create_dataset("features", data=np.asarray(features, dtype=h5py.string_dtype()))
                group.create_dataset("barcodes", data=np.asarray(cells, dtype=h5py.string_dtype()))
        staged.replace(output)
    return output
