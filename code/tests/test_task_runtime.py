"""Regression checks for task defaults, thin shell calls and validated batching."""
import contextlib
import importlib.util
import json
import os
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO / "code/functions"))
import runtime


def offline_fixture(root):
    producer = root / "producer.py"
    producer.write_text("# offline resource producer\n")
    annotation = "Chromosome\tStart\tEnd\tStrand\tGene\tTranscription_Start_Site\tTranscript_type\nchr1\t0\t10\t+\tG1\t1\tprotein_coding\n"
    (root / "genome_annotation.tsv").write_text(annotation)
    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    (root / "chromsizes.tsv").write_text("Chromosome\tStart\tEnd\n" + "".join(f"{name}\t0\t100\n" for name in chromosomes))
    (root / "source_manifest.tsv").write_text("file\tsha256\nfixture\tfixture\n")
    values = {"assembly": "GRCh38.p14", "annotation_release": "Ensembl 116",
              "chromosome_style": "UCSC", "primary_chromosomes": "25",
              "producer_script_sha256": runtime.sha256(producer), "annotation_rows": "1"}
    (root / "resource_validation.tsv").write_text("field\tvalue\n" + "".join(f"{key}\t{value}\n" for key, value in values.items()))
    with (root / "artifact_manifest.tsv").open("w") as handle:
        handle.write("file\tbytes\tsha256\n")
        for name in ("genome_annotation.tsv", "chromsizes.tsv", "source_manifest.tsv", "resource_validation.tsv"):
            path = root / name
            handle.write(f"{name}\t{path.stat().st_size}\t{runtime.sha256(path)}\n")
    return {"offline_resource_root": str(root), "offline_resource_producer": str(producer)}


class RuntimeTests(unittest.TestCase):
    def test_task_sources_have_no_argument_parsers_or_retired_forwarders(self):
        retired = ("fit_comparator_case_networks.R", "fit_comparator_case_networks.py",
                   "run_scode_worker.R", "log_message_path.R")
        for path in (REPO / "code").rglob("*"):
            if "tests" in path.parts or "vendor" in path.parts or path.suffix not in {".R", ".py"}:
                continue
            source = path.read_text()
            for token in ("commandArgs(", "run_task_cli", "run_function_cli", "parse_args(", *retired):
                self.assertNotIn(token, source, str(path))

    def test_retired_overrides_do_not_change_task_defaults(self):
        expected = runtime.rscript_executable("beeline")
        with patch.dict(os.environ, {"BEELINE_RSCRIPT": "/missing", "RSCRIPT_BIN": "/missing",
                                     "SCENICPLUS_CORES": "99", "OMP_NUM_THREADS": "99"}):
            self.assertEqual(runtime.rscript_executable("beeline"), expected)
            self.assertEqual(runtime.task_environment()["OMP_NUM_THREADS"], "1")

    def test_shell_entries_execute_each_step_once_in_order(self):
        wrappers = [REPO / f"0{number}_{name}.sh" for number, name in enumerate(
            ("beeline", "scmultisim", "gse274113", "gse162170", "gse194122", "plotting", "source_data"), 1)]
        wrappers += list((REPO / "code/requirements").glob("*.sh"))
        with tempfile.TemporaryDirectory() as directory:
            fixture = Path(directory) / "repo"
            fixture.mkdir()
            (fixture / "code/functions").mkdir(parents=True)
            (fixture / "code/functions/runtime.sh").write_text('run_task() { if [[ "$1" == *.R ]]; then Rscript --vanilla "$@"; else python3 "$1"; fi; }\n')
            fake = Path(directory) / "python3"
            fake.write_text(f"#!{sys.executable}\nimport json,os,sys\nprint(json.dumps([os.getcwd(),sys.argv[1:]]))\n")
            fake.chmod(0o755)
            fake_r = Path(directory) / "Rscript"
            fake_r.write_text(f"#!{sys.executable}\nimport json,os,sys\nprint(json.dumps([os.getcwd(),sys.argv[2:]]))\n")
            fake_r.chmod(0o755)
            env = dict(os.environ, PATH=directory + os.pathsep + os.environ["PATH"],
                       MULTICSN_PYTHON="/missing", BEELINE_PYTHON="/missing", FULL_MULTICSN_CORES="99")
            expected = {
                "01_beeline.sh": ["run_unified_benchmark.py", "promote_unified_batch.py", "run_core_scaling.R"],
                "02_scmultisim.sh": ["run_multiseed.R", "analyze_multiseed.R", "run_tf_gene_ablation.R"],
                "03_gse274113.sh": ["build_control_object.R", "prepare_control_candidates.R", "fit_control_timepoints.R", "fit_replicate_stability.R", "summarize_leave_one_replicate.R", "summarize_pairwise_fold_stability.R", "summarize_domain_decomposition.R", "summarize_domain_stability.R", "build_heldout_truth.R", "run_pando.R", "prepare_scenicplus_inputs.py", "run_pycistopic.py", "run_scenicplus.py", "export_endpoints.R", "summarize_network_catalog.R", "evaluate_tf_gene_comparators.R", "summarize_method_evidence.R"],
                "04_gse162170.sh": ["fit_network.R", "extract_cell_programs.R", "summarize_cell_programs.R", "analyze_key_gene_gsea.R", "run_gsea.R", "extract_locus_evidence.R", "materialize_locus_evidence.R", "screen_locus_candidates.R", "rank_loci_by_cell_type.R"],
                "05_gse194122.sh": ["fit_state_networks.R", "summarize_state_transitions.R", "extract_trajectory.R", "recompute_erythroid_gex_umap.R", "extract_donor_state_expression.R", "prepare_donor_deletion_manifests.R", "run_network_batch.R:run_donor_deletion", "summarize_donor_stable_rewiring.R", "analyze_state_regulomes.R", "prepare_cell_budget_manifests.R", "run_network_batch.R:run_cell_budget", "summarize_cell_budget_sensitivity.R"],
                "06_plotting.sh": ["fig2.R", "fig3.R", "figs7.R", "fig4.R", "figs8.R", "figs9.R", "fig5.R", "figs10.R", "figs1.R", "figs2.R", "figs3.R", "figs4.R", "figs5.R", "figs6.R", "figs11.R", "figs12.R", "verify_figures.py"],
                "07_source_data.sh": ["source_data.py"],
            }
            for wrapper in wrappers:
                with self.subTest(wrapper=wrapper.name):
                    relative = wrapper.relative_to(REPO)
                    target = fixture / relative
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(wrapper, target)
                    result = subprocess.check_output(["bash", str(target)], cwd=directory, env=env, text=True)
                    tasks = []
                    for line in result.splitlines():
                        cwd, args = json.loads(line)
                        self.assertEqual(Path(cwd).resolve(), fixture.resolve())
                        self.assertNotIn("-c", args)
                        self.assertNotIn("-e", args)
                        script = args[0]
                        self.assertTrue((REPO / script).is_file())
                        name = Path(script).name
                        if name == "run_network_batch.R":
                            name += ":" + args[1]
                        tasks.append(name)
                    self.assertEqual(tasks, expected.get(wrapper.name, [wrapper.stem + ".py"]))

    def test_shell_stops_before_later_steps_after_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "code/functions").mkdir(parents=True)
            (root / "code/functions/runtime.sh").write_text('run_task() { printf "%s\\n" "$1"; case "$1" in *analyze_multiseed.R) return 7;; esac; }\n')
            shutil.copy2(REPO / "02_scmultisim.sh", root / "02_scmultisim.sh")
            result = subprocess.run(["bash", str(root / "02_scmultisim.sh")],
                                    cwd=directory, text=True, capture_output=True)
            self.assertEqual(result.returncode, 7)
            self.assertNotIn("run_tf_gene_ablation.R", result.stdout)

    def test_python_function_process_keeps_sibling_imports_and_literal_settings(self):
        with tempfile.TemporaryDirectory(prefix="task's ") as directory:
            root = Path(directory)
            (root / "sibling.py").write_text("VALUE = 7\n")
            task = root / "task.py"
            task.write_text("import json\nfrom pathlib import Path\nfrom sibling import VALUE\n"
                            "def run(output, settings):\n"
                            "    Path(output).write_text(json.dumps([VALUE, settings]))\n"
                            "if __name__ == '__main__':\n    raise RuntimeError('task entry should not execute')\n")
            output = root / 'quoted "output".json'
            settings = {"nested": {"path": root / "input", "enabled": True, "optional": None}}
            command = runtime.python_function_command(task, "run", executable=sys.executable,
                                                      output=output, settings=settings)
            subprocess.run(command, check=True)
            self.assertEqual(json.loads(output.read_text()),
                             [7, {"nested": {"path": str(root / "input"), "enabled": True, "optional": None}}])

    def test_r_function_process_preserves_nested_and_quoted_settings(self):
        with tempfile.TemporaryDirectory(prefix="R task's ") as directory:
            root = Path(directory)
            task = root / "task.R"
            output = root / 'quoted "output".json'
            task.write_text('source("code/functions/runtime.R")\n'
                            'run <- function(output, settings) jsonlite::write_json(settings, output, auto_unbox=TRUE, null="null")\n'
                            'if (sys.nframe()==0L) stop("task entry should not execute")\n')
            settings = {"path": str(root / "输入"), "enabled": True, "optional": None, "workers": 3}
            command = runtime.r_function_command(task, "run", output=output, settings=settings)
            self.assertIn("-e", command)
            subprocess.run(command, cwd=REPO, check=True)
            self.assertEqual(json.loads(output.read_text()), settings)
            output.unlink()
            command = runtime.r_function_command(task, "unknown", output=output, settings=settings)
            result = subprocess.run(command, cwd=REPO, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(output.exists())

    def test_figure_validation_rejects_missing_and_stale_outputs(self):
        spec = importlib.util.spec_from_file_location("figure_check", REPO / "code/plotting/verify_figures.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            module.REPO = Path(directory)
            module.CHECKPOINT = Path(directory) / "start.json"
            module.CHECKPOINT.write_text(json.dumps({"started_ns": 2_000_000_000}))
            with self.assertRaisesRegex(RuntimeError, "Missing rendered panels"):
                module.verify_figures()
            outputs = [module.REPO / f"figures/{name}{panel}.pdf" for name, panels in module.PANELS.items() for panel in ("", *panels)]
            self.assertEqual(len(outputs), 64)
            panel_e = module.REPO / "figures/figs7e.pdf"
            self.assertIn(panel_e, outputs)
            outputs[0].parent.mkdir()
            for output in outputs:
                output.write_bytes(b"fixture")
                os.utime(output, ns=(1_000_000_000, 1_000_000_000))
            with self.assertRaisesRegex(RuntimeError, "panels not refreshed"):
                module.verify_figures()
            for output in outputs:
                os.utime(output, ns=(3_000_000_000, 3_000_000_000))
            panel_e.unlink()
            with self.assertRaisesRegex(RuntimeError, "figs7e.pdf"):
                module.verify_figures()
            panel_e.write_bytes(b"fixture")
            os.utime(panel_e, ns=(3_000_000_000, 3_000_000_000))
            module.verify_figures()
            self.assertFalse(module.CHECKPOINT.exists())

    def test_failed_child_stops_the_task(self):
        with self.assertRaises(subprocess.CalledProcessError):
            runtime.run_task([sys.executable, "-c", "raise SystemExit(7)"])

    def test_termination_stops_the_running_child(self):
        with tempfile.TemporaryDirectory() as directory:
            marker = Path(directory) / "child.pid"
            child = f"import os,time;from pathlib import Path;Path({str(marker)!r}).write_text(str(os.getpid()));time.sleep(30)"
            driver = f"import sys;sys.path.insert(0,{str(REPO / 'code/functions')!r});from runtime import run_task;run_task([sys.executable,'-c',{child!r}])"
            process = subprocess.Popen([sys.executable, "-c", driver], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            pid = None
            try:
                deadline = time.monotonic() + 5
                while not marker.is_file() and process.poll() is None and time.monotonic() < deadline:
                    time.sleep(0.02)
                self.assertTrue(marker.is_file())
                pid = int(marker.read_text())
                process.terminate()
                process.wait(timeout=5)
                with self.assertRaises(ProcessLookupError):
                    os.kill(pid, 0)
            finally:
                if process.poll() is None:
                    process.kill()
                    process.wait()
                if pid is not None:
                    try:
                        os.kill(pid, 9)
                    except ProcessLookupError:
                        pass

    def test_offline_resource_corruption_is_rejected(self):
        spec = importlib.util.spec_from_file_location("scenic_task", REPO / "code/gse274113/run_scenicplus.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            settings = offline_fixture(Path(directory))
            module.validate_offline_genome_resources(settings)
            (Path(directory) / "chromsizes.tsv").write_text("corrupt\n")
            with self.assertRaisesRegex(RuntimeError, "Stale offline genome artifact"):
                module.validate_offline_genome_resources(settings)

    def test_scenicplus_task_selection_preserves_timing_scope_and_defaults(self):
        spec = importlib.util.spec_from_file_location("scenic_profiles", REPO / "code/gse274113/run_scenicplus.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        original = {name: dict(settings) for name, settings in module.TASK_SETTINGS.items()}
        with patch.object(module, "run_scenicplus") as run:
            module.main()
            run.assert_called_once_with(**original["network"])
            run.reset_mock()
            module.run_timing()
            run.assert_called_once_with(**original["timing"])
        self.assertEqual(module.TASK_SETTINGS, original)

    def test_timing_benchmark_calls_the_actual_scenicplus_task(self):
        spec = importlib.util.spec_from_file_location("timing_driver", REPO / "code/gse274113/run_timing_benchmark.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            module.REPO = Path(directory)
            module.OUTPUT_ROOT = module.REPO / "timing"
            (module.REPO / "data/gse274113/objects/GSE274113_SCENICplus_full/cistopic_models").mkdir(parents=True)
            steps = list(module.scenicplus_steps(["day 7"]))
            self.assertEqual(len(steps), 1)
            self.assertEqual(steps[0][:2], ("all timepoints", "scope_A_downstream_eregulons"))
            self.assertEqual(steps[0][2][0], module.SCENICPLUS_PYTHON)
            self.assertEqual(steps[0][2][1], "-c")
            self.assertIn("run_scenicplus.py", steps[0][2][2])
            self.assertIn("run_timing", steps[0][2][2])
            fit = list(module.multi_csn_steps(["day 7"]))[0][2]
            self.assertEqual(fit[2], "-e")
            self.assertIn("code/functions/network_fit.R", fit[3])
            for literal in ('state="day 7"', "expected_cells=3304L", 'profile="gse274113"'):
                self.assertIn(literal, fit[3])



if __name__ == "__main__":
    unittest.main()
