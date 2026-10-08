"""Admission checks and rollback for replacing formal artifact directories."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

REPO = Path(__file__).absolute().parents[2]
SPEC = importlib.util.spec_from_file_location("promotion", REPO / "code/beeline/promote_unified_batch.py")
promotion = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(promotion)
import file_io


class PublicationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.stage, self.target = self.root / "stage", self.root / "formal"
        self.backup = self.root / "backup"
        for name in ("networks", "provenance", "runtime", "evaluation", "metrics"):
            (self.stage / name).mkdir(parents=True)
        for name in ("metrics", "core_scaling"):
            (self.target / name).mkdir(parents=True)
            (self.target / name / "old.tsv").write_text("original " + name)
        self.original = file_io.directory_hashes(self.target)
        self.scaling_inode = (self.target / "core_scaling/old.tsv").stat().st_ino
        methods = ("multiCSN", "PPCOR")
        rows = []
        for method in methods:
            for case in (1, 2):
                name = f"networks/{method}__case{case}.csv"
                (self.stage / name).write_text("regulator,target,weight\nA,B,1\n")
                rows.append(dict(Case=f"case{case}", Method=method, DatasetType="curated", DatasetGroup="HSC",
                                 RuntimeSeconds=float(case), TimingScope="inference_and_export", Status="complete",
                                 Network=name, NetworkSHA256=promotion.sha256(self.stage / name)))
        table = pd.DataFrame(rows)
        table.to_csv(self.stage / "provenance/run_manifest.tsv", sep="\t", index=False)
        table.drop(columns=["Network", "NetworkSHA256", "Status"]).to_csv(
            self.stage / "runtime/case_runtime.tsv", sep="\t", index=False)
        pd.DataFrame([dict(Method=method, cases=2, total_runtime_seconds=3.0, median_runtime_seconds=1.5,
                           maximum_runtime_seconds=2.0, runtime_ratio_to_multicsn=1.0) for method in methods]).to_csv(
            self.stage / "runtime/runtime_summary.tsv", sep="\t", index=False)
        (self.stage / "status.json").write_text(json.dumps(dict(state="validated", cases=4,
                                                               completed_methods=methods, smoke_test=False)))
        for name in ("environment.json", "input_code_sha256.json"):
            (self.stage / "provenance" / name).write_text('{"fixture": "recorded"}')
        evaluated = table[["Case", "Method", "DatasetType", "DatasetGroup"]].assign(
            AUPRC=0.5, EarlyPrecision=0.25, EPR=1.0)
        evaluated.to_csv(self.stage / "evaluation/case_metrics.tsv", sep="\t", index=False)
        table[["Case", "Method", "DatasetType"]].to_csv(
            self.stage / "evaluation/source_index.tsv", sep="\t", index=False)
        (self.stage / "evaluation/metric_policy.tsv").write_text("rule\tvalue\nversion\tfixture\n")
        for name in ("beeline_method_summary.tsv", "beeline_case_manifest.tsv"):
            evaluated.to_csv(self.stage / "evaluation" / name, sep="\t", index=False)
        for name in promotion.METRIC_FILES:
            evaluated.to_csv(self.stage / "metrics" / name,
                             sep="\t" if name.endswith(".tsv") else ",", index=False)
        ranks = evaluated.melt(id_vars=["Case", "Method", "DatasetType", "DatasetGroup"],
                               var_name="Metric", value_name="Value")
        ranks["Rank"] = 1.5
        ranks.to_csv(self.stage / "metrics/case_ranks.tsv", sep="\t", index=False)
        settings = patch.multiple(promotion, STAGE=self.stage, TARGET=self.target, BACKUP_ROOT=self.backup,
                                  METHODS=methods, EXPECTED_METHODS=2, EXPECTED_CASES=2, EXPECTED_RECORDS=4,
                                  case_paths=lambda limit: [Path(f"data/beeline/Curated/HSC/case{i}/ExpressionData.csv")
                                                            for i in (1, 2)])
        settings.start()
        self.addCleanup(settings.stop)

    def publish(self):
        with contextlib.redirect_stdout(io.StringIO()):
            promotion.main()

    def assert_original(self):
        self.assertEqual(file_io.directory_hashes(self.target), self.original)
        self.assertEqual((self.target / "core_scaling/old.tsv").stat().st_ino, self.scaling_inode)

    def test_success_preserves_scaling_and_verified_backup(self):
        self.publish()
        self.assertEqual((self.target / "core_scaling/old.tsv").stat().st_ino, self.scaling_inode)
        self.assertEqual(file_io.directory_hashes(self.target / "metrics"),
                         file_io.directory_hashes(self.stage / "metrics"))
        self.assertEqual(file_io.directory_hashes(self.target / "evaluation"),
                         file_io.directory_hashes(self.stage / "evaluation"))
        self.assertEqual(len(list((self.target / "networks").glob("*.csv"))), 4)
        backup = next(self.backup.iterdir())
        self.assertEqual(file_io.directory_hashes(backup / "contents"), json.loads((backup / "sha256.json").read_text()))
        self.assertEqual((backup / "contents/metrics/old.tsv").read_text(), "original metrics")

    def test_missing_required_file_leaves_original_intact(self):
        (self.stage / "runtime/runtime_summary.tsv").unlink()
        with self.assertRaisesRegex(SystemExit, "missing required"):
            self.publish()
        self.assert_original()

    def test_invalid_runtime_case_time_and_summary_are_rejected(self):
        path = self.stage / "runtime/case_runtime.tsv"
        original = path.read_bytes()
        for column, value in (("Case", "unrelated"), ("RuntimeSeconds", -1), ("RuntimeSeconds", 1.1)):
            with self.subTest(column=column, value=value):
                path.write_bytes(original)
                table = pd.read_csv(path, sep="\t")
                table.loc[0, column] = value
                table.to_csv(path, sep="\t", index=False)
                with self.assertRaises(SystemExit):
                    self.publish()
                self.assert_original()
        path.write_bytes(original)
        summary = self.stage / "runtime/runtime_summary.tsv"
        summary.write_text(summary.read_text().replace("3.0", "4.0"))
        with self.assertRaisesRegex(SystemExit, "summary disagrees"):
            self.publish()
        self.assert_original()

    def test_missing_or_inconsistent_evaluation_cannot_be_promoted(self):
        path = self.stage / "evaluation/case_metrics.tsv"
        original = path.read_bytes()
        for column, value in (("Case", "unregistered"), ("EarlyPrecision", 0.75)):
            with self.subTest(column=column):
                table = pd.read_csv(path, sep="\t")
                table.loc[0, column] = value
                table.to_csv(path, sep="\t", index=False)
                with self.assertRaisesRegex(SystemExit, "tables disagree"):
                    self.publish()
                self.assert_original()
                path.write_bytes(original)
        path.unlink()
        with self.assertRaisesRegex(SystemExit, "missing required"):
            self.publish()
        self.assert_original()

    def test_copy_failure_leaves_original_intact(self):
        copy = promotion.shutil.copy2
        def fail_copy(source, destination, *args, **kwargs):
            if Path(source).name == "multiCSN__case2.csv":
                raise OSError("injected copy failure")
            return copy(source, destination, *args, **kwargs)
        with patch.object(promotion.shutil, "copy2", side_effect=fail_copy):
            with self.assertRaisesRegex(OSError, "injected"):
                self.publish()
        self.assert_original()

    def test_unregistered_case_inventory_is_rejected(self):
        with patch.object(promotion, "case_paths", return_value=[
                Path(f"data/beeline/Curated/HSC/other{i}/ExpressionData.csv") for i in (1, 2)]):
            with self.assertRaisesRegex(SystemExit, "case coverage"):
                self.publish()
        self.assert_original()

    def test_cleanup_failure_reports_successful_publication(self):
        remove = file_io.shutil.rmtree
        def fail_cleanup(path, *args, **kwargs):
            if Path(path).name.startswith(".formal-rollback-"):
                raise PermissionError("injected cleanup failure")
            return remove(path, *args, **kwargs)
        with patch.object(file_io.shutil, "rmtree", side_effect=fail_cleanup), contextlib.redirect_stderr(io.StringIO()) as log:
            self.publish()
        self.assertIn("old temporary tree retained", log.getvalue())
        self.assertEqual(len(list((self.target / "networks").glob("*.csv"))), 4)
        self.assertEqual((self.target / "core_scaling/old.tsv").stat().st_ino, self.scaling_inode)

    def test_switch_failure_rolls_back_original_and_scaling(self):
        rename = Path.rename
        def fail_switch(path, destination):
            if path.name == "results" and destination == self.target:
                raise OSError("injected switch failure")
            return rename(path, destination)
        with patch.object(Path, "rename", fail_switch):
            with self.assertRaisesRegex(OSError, "injected"):
                self.publish()
        self.assert_original()

    def test_corrupted_backup_prevents_switch(self):
        copy = file_io.shutil.copy2
        def corrupt_backup(source, destination, *args, **kwargs):
            result = copy(source, destination, *args, **kwargs)
            if "contents" in Path(destination).parts:
                Path(destination).write_text("injected corruption")
            return result
        # copytree uses its bound copy_function, so corrupt its explicit argument.
        tree = file_io.shutil.copytree
        def copy_tree(source, destination, **kwargs):
            return tree(source, destination, copy_function=corrupt_backup, **kwargs)
        with patch.object(file_io.shutil, "copytree", side_effect=copy_tree):
            with self.assertRaisesRegex(RuntimeError, "Backup verification failed"):
                self.publish()
        self.assert_original()


if __name__ == "__main__":
    unittest.main()
