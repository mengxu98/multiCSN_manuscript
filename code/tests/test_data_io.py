"""Contracts for reusable file manifests and sparse HDF5 input readers."""
import csv
import io
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "functions"))
from file_io import write_artifact_manifest
from hdf5_io import decode_hdf5_strings, read_sparse_group
from runtime import sha256


class DataIOTests(unittest.TestCase):
    def test_hash_and_manifest_coverage(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "nested").mkdir()
            (root / "a.tsv").write_bytes(b"abc")
            (root / "nested" / 'gene "A".tsv').write_bytes(b"")
            write_artifact_manifest(root)
            first = (root / "artifact_manifest.tsv").read_bytes()
            write_artifact_manifest(root)
            self.assertEqual((root / "artifact_manifest.tsv").read_bytes(), first)
            rows = list(csv.DictReader(io.StringIO(first.decode()), delimiter="\t"))
            self.assertEqual([row["file"] for row in rows], ["a.tsv", 'nested/gene "A".tsv'])
            self.assertEqual([int(row["bytes"]) for row in rows], [3, 0])
            self.assertEqual(sha256(root / "a.tsv"),
                             "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
            for row in rows:
                self.assertEqual(row["sha256"], sha256(root / row["file"]))

    def test_decoded_labels(self):
        self.assertEqual(decode_hdf5_strings([b"cell1", "cell2", 3]), ["cell1", "cell2", "3"])

    def test_sparse_orientation_and_empty_group(self):
        with h5py.File("input.h5", "w", driver="core", backing_store=False) as handle:
            for name, data, indices, indptr in (
                ("counts", [1.0, 3.0, 2.0], [0, 1, 0], [0, 1, 2, 3]),
                ("empty", [], [], [0, 0, 0, 0]),
            ):
                group = handle.create_group(name)
                for key, values in (("data", data), ("indices", indices), ("indptr", indptr),
                                    ("shape", [2, 3]), ("features", [b"f1", b"f2"]),
                                    ("barcodes", [b"c1", b"c2", b"c3"])):
                    group.create_dataset(key, data=values)
                matrix, features, barcodes = read_sparse_group(handle, name)
                self.assertEqual(matrix.shape, (2, 3))
                self.assertEqual(features, ["f1", "f2"])
                self.assertEqual(barcodes, ["c1", "c2", "c3"])
                expected = [[1, 0, 2], [0, 3, 0]] if name == "counts" else [[0, 0, 0], [0, 0, 0]]
                np.testing.assert_array_equal(matrix.toarray(), expected)


if __name__ == "__main__":
    unittest.main()
