"""Interpreter locations and process settings used by analysis tasks."""
import contextlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time

sys.dont_write_bytecode = True

REPO = Path(__file__).absolute().parents[2]
EXECUTION_ROOT = REPO.parent
USER_HOME = Path.home()
PYTHON_LOCATIONS = {
    "general": (USER_HOME / ".omicos/env/.venv/bin/python",),
    "beeline": (USER_HOME / "miniconda3/envs/beeline_py39/bin/python",
                USER_HOME / ".omicos/env/.venv/bin/python"),
    "grnboost2": (USER_HOME / ".venvs/beeline-unified-py39/bin/python",
                  USER_HOME / ".omicos/env/grnboost2-py39/bin/python"),
    "scenicplus": (EXECUTION_ROOT / "comparator_envs/scenicplus/bin/python",),
    "linger": (EXECUTION_ROOT / "comparator_envs/linger/bin/python",),
}
R_LOCATIONS = {
    "analysis": (USER_HOME / "miniconda3/envs/R_4.5/bin/Rscript",),
    "beeline": (USER_HOME / "miniconda3/envs/R_4.5/bin/Rscript",),
    "plotting": (USER_HOME / "miniconda3/envs/R_4.5/bin/Rscript",),
}
THREAD_SETTINGS = {key: "1" for key in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLAS_NUM_THREADS",
    "BLIS_NUM_THREADS",
)}


def python_executable(role="general"):
    for candidate in PYTHON_LOCATIONS[role]:
        if candidate.is_file():
            return str(candidate)
    if role in {"general", "beeline"}:
        return sys.executable
    raise RuntimeError(f"Missing {role} Python environment; see code/requirements/install_environments.py")


def rscript_executable(role="analysis"):
    for candidate in R_LOCATIONS[role]:
        if candidate.is_file():
            return str(candidate)
    candidate = shutil.which("Rscript")
    if candidate:
        return candidate
    raise RuntimeError("Rscript is unavailable; see code/requirements/install_environments.py")


def r_function_command(script, function, *, executable=None, role="analysis", **settings):
    def literal(value):
        if value is None:
            return "NULL"
        if isinstance(value, bool):
            return "TRUE" if value else "FALSE"
        if isinstance(value, int):
            return f"{value}L"
        if isinstance(value, dict):
            return "list(" + ", ".join(f"{key}={literal(item)}" for key, item in value.items()) + ")"
        if isinstance(value, (list, tuple)):
            return "list(" + ", ".join(literal(item) for item in value) + ")"
        return json.dumps(str(value) if isinstance(value, Path) else value)
    arguments = ", ".join(f"{key}={literal(value)}" for key, value in settings.items())
    return [executable or rscript_executable(role), "--vanilla", "-e",
            f"source({literal(REPO / script)}); {function}({arguments})"]


def python_function_command(script, function, *, executable=None, **settings):
    path = REPO / script
    def plain(value):
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, dict):
            return {key: plain(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [plain(item) for item in value]
        return value
    expression = (f"import runpy,sys; sys.path.insert(0, {str(path.parent)!r}); "
                  f"module=runpy.run_path({str(path)!r}); module[{function!r}](**{plain(settings)!r})")
    return [executable or python_executable(), "-B", "-c", expression]


def task_environment(python=None, r_library=None, hash_seed=None, mallet_memory=None):
    env = dict(os.environ)
    env.update(THREAD_SETTINGS, LC_ALL="C", PYTHONDONTWRITEBYTECODE="1")
    if python:
        env["PATH"] = str(Path(python).parent) + os.pathsep + env.get("PATH", "")
    if r_library is not None:
        env["R_LIBS_USER"] = str(r_library)
    if hash_seed is not None:
        env["PYTHONHASHSEED"] = str(hash_seed)
    if mallet_memory is not None:
        env["MALLET_MEMORY"] = str(mallet_memory)
    return env


def use_python(role="general", hash_seed=None, mallet_memory=None):
    """Select the declared environment before a task imports numerical packages."""
    executable = python_executable(role)
    env = task_environment(executable, hash_seed=hash_seed, mallet_memory=mallet_memory)
    if Path(sys.executable).absolute() != Path(executable).absolute() or (
        hash_seed is not None and os.environ.get("PYTHONHASHSEED") != str(hash_seed)
    ):
        os.execve(executable, [executable, *sys.argv], env)
    os.environ.update(env)


def utc_now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def raise_on_termination(signum, frame):
    raise KeyboardInterrupt


@contextlib.contextmanager
def termination_guard():
    previous = None
    if threading.current_thread() is threading.main_thread():
        previous = signal.signal(signal.SIGTERM, raise_on_termination)
    try:
        yield
    finally:
        if previous is not None:
            signal.signal(signal.SIGTERM, previous)


def run_task(command, env=None, cwd=REPO, stage=None, workers=1,
             monitor=False, timed=False):
    """Run a task with optional Linux timing and process-tree RSS sampling."""
    env = task_environment() if env is None else env
    command = [str(value) for value in command]
    if timed or monitor:
        command = ["/usr/bin/time", "-v", *command]
    with termination_guard(), tempfile.TemporaryDirectory(prefix="multicsn-resource-") as directory:
        output = Path(directory) / "rss.tsv"
        sampler = None
        if monitor:
            print(f"resource_stage\t{stage}\nresource_requested_workers\t{workers}", flush=True)
        process = subprocess.Popen(command, cwd=cwd, env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, start_new_session=True)
        try:
            if monitor:
                sampler = subprocess.Popen(python_function_command(
                    "code/functions/sample_process_tree_rss.py", "sample_rss", executable=sys.executable,
                    root_pid=process.pid, output=output, interval=0.25), env=env)
            for line in process.stdout:
                print(line, end="", flush=True)
            status = process.wait()
            sampler_status = sampler.wait() if sampler else 0
        except BaseException:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                process.wait()
            if sampler and sampler.poll() is None:
                sampler.terminate()
                sampler.wait()
            raise
        finally:
            process.stdout.close()
        if monitor:
            if output.is_file() and output.stat().st_size:
                print(output.read_text(), end="", flush=True)
            else:
                print("process_tree_rss_sum_peak_kb\tNA\nprocess_tree_processes_at_peak\tNA\n"
                      "process_tree_rss_samples\t0\nprocess_tree_rss_interval_seconds\t0.25")
            print(f"resource_command_exit_status\t{status}\nresource_sampler_exit_status\t{sampler_status}", flush=True)
        if status or sampler_status:
            raise subprocess.CalledProcessError(status or sampler_status, command)
