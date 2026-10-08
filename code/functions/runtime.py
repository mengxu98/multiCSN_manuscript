"""Interpreter locations and process settings used by analysis tasks."""
import contextlib
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import threading
import time

REPO = Path(__file__).absolute().parents[2]
EXECUTION_ROOT = REPO.parent
USER_HOME = Path.home()
PYTHON_LOCATIONS = {
    "general": (USER_HOME / ".omicos/env/.venv/bin/python",),
    "beeline": (USER_HOME / "miniconda3/envs/beeline_py39/bin/python",),
    "grnboost2": (USER_HOME / ".venvs/beeline-unified-py39/bin/python",
                  USER_HOME / ".omicos/env/grnboost2-py39/bin/python"),
    "scenicplus": (EXECUTION_ROOT / "comparator_envs/scenicplus/bin/python",),
    "linger": (EXECUTION_ROOT / "comparator_envs/linger/bin/python",),
}
R_LOCATIONS = {
    "analysis": (USER_HOME / "miniconda3/envs/R_4.4/bin/Rscript",),
    "beeline": (USER_HOME / "miniconda3/envs/beeline_r453/bin/Rscript",),
    "plotting": (USER_HOME / "miniconda3/envs/R_4.4/bin/Rscript",),
}
COMPARATOR_R_LIBRARY = EXECUTION_ROOT / "comparator_r_library"
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
    return [executable or python_executable(), "-c", expression]


def task_environment(python=None, r_library=None, hash_seed=None, mallet_memory=None):
    env = dict(os.environ)
    env.update(THREAD_SETTINGS, LC_ALL="C")
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


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


class _Tee:
    def __init__(self, original, handle):
        self.original, self.handle = original, handle

    def write(self, value):
        self.original.write(value)
        self.handle.write(value)
        self.flush()

    def flush(self):
        self.original.flush()
        self.handle.flush()


@contextlib.contextmanager
def tee_log(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        stdout, stderr = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = _Tee(stdout, handle), _Tee(stderr, handle)
        try:
            yield
        finally:
            sys.stdout, sys.stderr = stdout, stderr


def run_task(command, env=None, cwd=REPO, stage=None, workers=1, producer=None,
             monitor=False, timed=False):
    """Run a task with optional Linux timing and process-tree RSS sampling."""
    env = task_environment() if env is None else env
    command = [str(value) for value in command]
    if monitor and (not producer or not Path(producer).is_file() or workers < 1):
        raise ValueError("Invalid resource producer or worker count")
    if timed or monitor:
        command = ["/usr/bin/time", "-v", *command]
    with termination_guard(), tempfile.TemporaryDirectory(prefix="multicsn-resource-") as directory:
        output = Path(directory) / "rss.tsv"
        sampler = None
        if monitor:
            machine_hash = sha256("/etc/machine-id")
            model = next((line.split(":", 1)[1].strip() for line in
                          Path("/proc/cpuinfo").read_text().splitlines()
                          if line.startswith("model name")), "")
            host = socket.getfqdn()
            kernel = os.uname().release
            cpus = os.cpu_count()
            fingerprint = hashlib.sha256(
                f"{host}\t{machine_hash}\t{cpus}\t{model}\t{kernel}\n".encode()).hexdigest()
            run_id = f"{stage}_{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}_{os.getpid()}"
            fields = {
                "resource_block_started": run_id, "resource_block_started_utc": utc_now(),
                "resource_stage": stage, "resource_requested_workers": workers,
                "resource_host_name": host, "resource_machine_id_sha256": machine_hash,
                "resource_boot_id": Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
                "resource_host_fingerprint": fingerprint, "resource_cpu_online": cpus,
                "resource_cpu_model": model, "resource_kernel_release": kernel,
                "resource_producer_path": producer, "resource_producer_sha256": sha256(producer),
                "resource_wrapper_sha256": sha256(__file__),
                "resource_sampler_sha256": sha256(REPO / "code/functions/sample_process_tree_rss.py"),
                "resource_thread_environment": ";".join(f"{key}={env[key]}" for key in THREAD_SETTINGS),
                "resource_command": shlex.join(command[2:]),
                "resource_command_started_epoch_ns": time.time_ns(),
            }
            for key, value in fields.items():
                print(f"{key}\t{value}", flush=True)
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
            print(f"resource_command_exit_status\t{status}\nresource_sampler_exit_status\t{sampler_status}\n"
                  f"resource_command_completed_epoch_ns\t{time.time_ns()}\n"
                  f"resource_block_completed_utc\t{utc_now()}\nresource_block_completed\t{run_id}", flush=True)
        if status or sampler_status:
            raise subprocess.CalledProcessError(status or sampler_status, command)
