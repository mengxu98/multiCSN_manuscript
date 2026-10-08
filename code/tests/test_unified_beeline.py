"""Contract tests; do not invoke scientific inference or touch formal data."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

SPEC = importlib.util.spec_from_file_location(
    "unified", Path(__file__).resolve().parents[1] / "beeline/run_unified_benchmark.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class UnifiedContractTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for name in ("work", "provenance", "networks"):
            (self.root / name).mkdir(parents=True)
        self.net = self.root / "networks/multiCSN__case.csv"
        self.net.write_text("regulator,target,weight\na,b,1\n")
        self.row = dict(Case="case", Method="multiCSN", DatasetType="curated", DatasetGroup="HSC",
                        RuntimeSeconds=0.1, Network=str(self.net.relative_to(self.root)),
                        NetworkSHA256=MODULE.sha256(self.net), TimingScope="inference_and_export",
                        Status="complete")
        self.write_rows([self.row])
        self.paths = patch.object(MODULE, "case_paths", return_value=[
            Path("data/beeline/Curated/HSC/case/ExpressionData.csv")])
        self.paths.start()
        self.addCleanup(self.paths.stop)

    def write_rows(self, rows):
        pd.DataFrame(rows).to_csv(self.root / "work/multiCSN.tsv", sep="\t", index=False)

    def validate(self):
        return MODULE.validate(self.root, ["multiCSN"], 1)

    def test_paired_network_runtime(self):
        self.assertEqual(self.validate(), 1)
        self.assertEqual(sorted(p.name for p in (self.root / "runtime").iterdir()),
                         ["case_runtime.tsv", "runtime_summary.tsv"])

    def test_changed_network_rejected(self):
        self.net.write_text("regulator,target,weight\na,b,2\n")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            self.validate()

    def test_duplicate_case_rejected(self):
        self.write_rows([self.row, self.row])
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            self.validate()

    def test_zero_time_rejected(self):
        self.row["RuntimeSeconds"] = 0
        self.write_rows([self.row])
        with self.assertRaisesRegex(ValueError, "Invalid timings"):
            self.validate()

    def test_case_network_mismatch_rejected(self):
        self.row["Network"] = "networks/multiCSN__other.csv"
        self.write_rows([self.row])
        with self.assertRaisesRegex(ValueError, "identity mismatch"):
            self.validate()

    def test_unfinished_record_rejected(self):
        self.row["Status"] = "running"
        self.write_rows([self.row])
        with self.assertRaisesRegex(ValueError, "Unfinished"):
            self.validate()

    def test_input_order_matches_r_full_path_sort(self):
        # Path comparison sorts components, unlike R sort(full_path).
        self.paths.stop()
        inputs = [Path(f"/repo/data/beeline/Curated/GSD/GSD-2000-{n}/ExpressionData.csv")
                  for n in range(1, 501)]
        with patch.object(Path, "rglob", side_effect=[inputs, []]), patch.object(MODULE, "REPO", Path("/repo")):
            self.assertEqual(MODULE.case_paths(1), sorted(inputs, key=str)[:1])


if __name__ == "__main__":
    unittest.main()
