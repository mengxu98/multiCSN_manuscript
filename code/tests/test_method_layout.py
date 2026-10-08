"""Check method imports and repository paths without running inference."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import pandas as pd

REPO = Path(__file__).resolve().parents[2]


def load_driver(name):
    spec = importlib.util.spec_from_file_location(name, REPO / "code/beeline" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MethodLayoutTests(unittest.TestCase):
    def test_drivers_resolve_the_repository(self):
        for name in ("run_unified_benchmark", "run_grnboost2"):
            with self.subTest(driver=name):
                self.assertEqual(load_driver(name).REPO, REPO)

    def test_python_dispatch_uses_method_modules(self):
        load_driver("run_grnboost2")
        from methods.GNIPLR import infer_gniplr
        from methods.SCSGL import infer_scsgl
        for name, runner in {"GNIPLR": infer_gniplr, "SCSGL": infer_scsgl}.items():
            self.assertEqual(runner.__module__, f"methods.{name}")
        from methods.SCSGL import REPO as scsgl_repo
        self.assertTrue((scsgl_repo / "code/vendor/SCSGL/run_scSGL.py").is_file())

    def test_grnboost2_candidate_selection(self):
        load_driver("run_grnboost2")
        from methods.GRNBoost2 import infer_tf_names
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ExpressionData.csv"
            expression = pd.DataFrame(columns=["a", "b", "c"])
            self.assertEqual(infer_tf_names(expression, path), ["a", "b", "c"])
            pd.DataFrame({"RequestedTFs": [2]}).to_csv(path.with_name("case_metadata.csv"), index=False)
            self.assertEqual(infer_tf_names(expression, path), ["a", "b"])
            expression.columns = ["TF1", "b", "TF2"]
            self.assertEqual(infer_tf_names(expression, path), ["TF1", "TF2"])


if __name__ == "__main__":
    unittest.main()
