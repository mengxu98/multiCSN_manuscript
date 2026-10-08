"""Check the actual Source Data loader against the plotted locus identities."""
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("source_data", ROOT / "code/source_data/source_data.py")
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)


class LocusSources(unittest.TestCase):
    def test_stability_source_labels_follow_current_figure_panels(self):
        expected = {
            "FigS7_fold_stability": "Figure S7c",
            "FigS7_pairwise_stability": "Figure S7c",
            "FigS7_top100": "Figure S7d",
            "FigS7_rank_deciles": "Figure S7e",
            "FigS7e_endpoint_counts": "Figure S7b",
        }
        entries = {entry["sheet"]: entry for entry in source.SOURCE_SHEETS}
        for sheet, panel in expected.items():
            self.assertEqual(entries[sheet]["figure"], panel)
            self.assertTrue((ROOT / entries[sheet]["path"]).is_file())

    def test_path_identity_and_count(self):
        for gene, count in (("ZNF804A", 32), ("KCNC2", 33)):
            entry = next(e for e in source.SOURCE_SHEETS if gene in e["sheet"] and e["sheet"].endswith("_paths"))
            table = source.load_source_sheet(entry)
            self.assertEqual(set(table["Target"]), {gene})
            self.assertEqual(len(table), count)

    def test_wrong_target_rejected(self):
        entry = next(e for e in source.SOURCE_SHEETS if "ZNF804A" in e["sheet"] and e["sheet"].endswith("_paths"))
        wrong_target = source.pd.DataFrame({"target": ["ADRA2A"] * 32})
        with patch.object(source, "read_table", return_value=wrong_target):
            with self.assertRaises(ValueError):
                source.load_source_sheet(entry)

    def test_all_locus_tables_use_current_candidates(self):
        entries = [e for e in source.SOURCE_SHEETS if any(g in e["sheet"] for g in ("ZNF804A", "KCNC2"))]
        self.assertEqual(len(entries), 8)
        for entry in entries:
            self.assertIn("/loci/candidates/", entry["path"])
            self.assertFalse(source.load_source_sheet(entry).empty)
            if "KCNC2" in entry["sheet"]:
                self.assertNotIn("Figure 4", entry["figure"])


if __name__ == "__main__":
    unittest.main()
