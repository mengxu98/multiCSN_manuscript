#!/usr/bin/env python3
"""Install the declared analysis environments; The dry_run setting prints commands only."""
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "functions"))
from runtime import REPO, EXECUTION_ROOT, USER_HOME

CONDA = USER_HOME / "miniconda3/bin/conda"
R_ENV_NAME = "beeline_r453"
R_ENV = USER_HOME / "miniconda3/envs" / R_ENV_NAME
PY_ENV = USER_HOME / ".omicos/env/.venv"
ENV_ROOT = EXECUTION_ROOT / "comparator_envs"
SOURCE_ROOT = EXECUTION_ROOT / "sources"
WORK = USER_HOME / "srcpkgs"
GRNBOOST2_ENV = USER_HOME / ".venvs/beeline-unified-py39"
PYPI = "https://pypi.org/simple"
CLONE_RETRIES = 4


class Installer:
    def __init__(self, dry_run=False, force=False):
        self.dry_run, self.force = dry_run, force

    def run(self, command, **kwargs):
        command = [str(value) for value in command]
        if self.dry_run:
            print("[dry-run] " + shlex.join(command))
        else:
            subprocess.run(command, cwd=REPO, check=True, **kwargs)

    def clone(self, url, destination):
        if (destination / ".git").is_dir() and not self.force:
            return
        if self.dry_run:
            self.run(["git", "clone", "-q", "--depth", "1", url, destination])
            return
        destination.parent.mkdir(parents=True, exist_ok=True)
        for attempt in range(CLONE_RETRIES):
            shutil.rmtree(destination, ignore_errors=True)
            result = subprocess.run(["git", "clone", "-q", "--depth", "1", url, str(destination)])
            if result.returncode == 0:
                return
            if attempt + 1 < CLONE_RETRIES:
                time.sleep(10)
        raise RuntimeError(f"Unable to clone {url}")

    def available(self, python, modules):
        return not self.force and Path(python).is_file() and subprocess.run(
            [str(python), "-c", modules], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0

    def r(self):
        if not R_ENV.is_dir() or self.force:
            self.run([CONDA, "create", "-y", "-n", R_ENV_NAME, "-c", "conda-forge", "r-base=4.5", "r-essentials"])
        rscript = R_ENV / "bin/Rscript"
        expression = "pkgs <- readLines(\"code/requirements/r-requirements.txt\"); pkgs <- pkgs[nzchar(pkgs) & !grepl('^#', pkgs)]; pkgs <- pkgs[!grepl('/', pkgs)]; need <- pkgs[!vapply(pkgs, function(p) requireNamespace(p, quietly=TRUE), logical(1))]; if (length(need)) install.packages(need, repos='https://cloud.r-project.org', Ncpus=4)"
        self.run([rscript, "--vanilla", "-e", expression])
        self.run([rscript, "--vanilla", REPO / "code/requirements/install_csn_packages.R"])
        self.run([rscript, "--vanilla", "-e", "remotes::install_github(c('mengxu98/thisutils','mengxu98/thisplot'), dependencies=NA, upgrade='never')"])
        self.run([rscript, "--vanilla", "-e", "pkgs <- c('inferCSN','multiCSN','GENIE3','LEAP','SCORPION','ppcor','glmnet','kSamples','thisutils','thisplot'); ok <- vapply(pkgs, function(p) requireNamespace(p, quietly=TRUE), logical(1)); if (any(!ok)) stop('Missing packages: ', paste(pkgs[!ok], collapse=',')); stopifnot('greedy_l0' %in% eval(formals(multiCSN::fit_model)$method))"])

    def py(self):
        if not (PY_ENV / "bin/python").is_file():
            self.run([sys.executable, "-m", "venv", PY_ENV])
        python = PY_ENV / "bin/python"
        self.run([python, "-m", "pip", "install", "-q", "--index-url", PYPI, "--upgrade", "pip"])
        self.run([python, "-m", "pip", "install", "-q", "--index-url", PYPI, "-r", REPO / "code/requirements/py-requirements.txt"])
        self.run([python, "-c", "import pandas, openpyxl; print('pandas/openpyxl OK')"])

    def scenicplus(self):
        prefix = ENV_ROOT / "scenicplus"
        python = prefix / "bin/python"
        if self.available(python, "import pycisTopic, pycistarget, scenicplus"):
            print(f"Already available: {prefix}")
            return
        source = SOURCE_ROOT / "scenicplus"
        if not (source / ".git").is_dir():
            self.clone("https://github.com/aertslab/SCENICplus.git", source)
            if not self.dry_run:
                commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True)
                record = REPO / "logs/gse274113/comparator_install/scenicplus_source_commit.txt"
                record.parent.mkdir(parents=True, exist_ok=True)
                record.write_text(commit)
        packages = "python=3.11.8 pip numpy scipy scikit-learn pandas pyarrow numba cython pytables h5py pysam bedtools macs2 gensim datrie statsmodels patsy seaborn annoy python-igraph leidenalg louvain fa2 nmslib meson-python meson ninja pkg-config scikit-build-core cmake".split()
        self.run([CONDA, "create", "-y", "-p", prefix, "-c", "conda-forge", "-c", "bioconda", *packages])
        self.run([python, "-m", "pip", "install", "-q", "--index-url", PYPI, "setuptools<81"])
        for name in ("LoomXpy", "pycisTopic", "pycistarget"):
            target = WORK / name
            self.clone(f"https://github.com/aertslab/{name}.git", target)
            flags = [] if name == "LoomXpy" else ["--no-build-isolation"]
            self.run([python, "-m", "pip", "install", "-q", "--index-url", PYPI, *flags, target])
        self.run([python, "-m", "pip", "install", "-q", "--index-url", PYPI, "--no-deps", "--no-build-isolation", source])
        self.run([python, "-m", "pip", "install", "-q", "--index-url", PYPI, *"snakemake polars pyranges scatac_fragment_tools ray matplotlib adjustText lda tmtoolkit".split()])
        self.run([python, "-c", "from pycisTopic.cistopic_class import create_cistopic_object; from pycisTopic.lda_models import run_cgs_models_mallet; from pycisTopic.topic_binarization import binarize_topics; from pycisTopic.utils import region_names_to_coordinates; import pycistarget, scenicplus, anndata, scanpy, snakemake, pandas, numpy"])
        self.run([python, "-m", "pip", "freeze"])

    def linger(self):
        prefix = ENV_ROOT / "linger"
        python = prefix / "bin/python"
        if not python.is_file():
            self.run([CONDA, "create", "-y", "-p", prefix, "-c", "conda-forge", "-c", "bioconda", "python=3.10", "pip", "pysam", "bedtools", "h5py"])
        version_check = "from importlib.metadata import version; assert version('LingerGRN') == '1.110'"
        if not self.available(python, version_check):
            flags = ["--force-reinstall"] if self.force else []
            self.run([python, "-m", "pip", "install", "-q", "--index-url", PYPI, *flags, "LingerGRN==1.110", "psutil"])
        self.clone("https://github.com/Durenlab/LINGER.git", SOURCE_ROOT / "LINGER")
        if self.dry_run:
            print("[dry-run] Apply LingerGRN-1.110-trans_shap-missing-key.patch to LingerGRN/LL_net.py")
        else:
            locate = (version_check + "; from importlib.util import find_spec; from pathlib import Path; "
                      "print(Path(find_spec('LingerGRN.LL_net').origin).parent)")
            site = Path(subprocess.check_output([str(python), "-c", locate], text=True).strip())
            patch = REPO / "code/gse274113/LingerGRN-1.110-trans_shap-missing-key.patch"
            if not (site / "LL_net.py").is_file():
                raise RuntimeError("LingerGRN/LL_net.py is missing")
            if subprocess.run(["patch", "--batch", "--forward", "-p0", "--dry-run", "-i", str(patch)], cwd=site, capture_output=True).returncode == 0:
                subprocess.run(["patch", "--batch", "--forward", "-p0", "-i", str(patch)], cwd=site, check=True)
                for cache in site.rglob("__pycache__"):
                    shutil.rmtree(cache)
            elif subprocess.run(["patch", "--batch", "-R", "-p0", "--dry-run", "-i", str(patch)], cwd=site, capture_output=True).returncode:
                raise RuntimeError("LINGER patch does not match the installed source")
        self.run([python, "-c", "import LingerGRN, LingerGRN.LL_net, LingerGRN.LINGER_tr, LingerGRN.preprocess, LingerGRN.pseudo_bulk"])

    def grnboost2(self):
        # The lock targets Python 3.9; the bootstrap interpreter may be newer.
        if not (GRNBOOST2_ENV / "bin/python").is_file() or self.force:
            self.run([CONDA, "create", "-y", "-p", GRNBOOST2_ENV, "-c", "conda-forge", "python=3.9", "pip"])
        python = GRNBOOST2_ENV / "bin/python"
        self.run([python, "-m", "pip", "install", "-q", "--upgrade", "pip"])
        self.run([python, "-m", "pip", "install", "-q", "--index-url", PYPI, "-r", REPO / "code/requirements/grnboost2-official-py39.lock.txt"])


SETTINGS = dict(parts=("all",), dry_run=False, force=False)


def main():
    selected = SETTINGS["parts"]
    names = ("r", "py", "scenicplus", "linger", "grnboost2")
    if set(selected) - set(names) - {"all"}:
        raise ValueError("Unknown environment selection")
    installer = Installer(SETTINGS["dry_run"], SETTINGS["force"])
    for name in names:
        if "all" in selected or name in selected:
            print(f"Installing {name}", flush=True)
            getattr(installer, name)()

if __name__ == "__main__":
    main()
