"""Exercise resource preparation and LINGER repair without network or real environments."""
import contextlib
import gzip
import importlib.util
import io
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import venv

REPO = Path(__file__).absolute().parents[2]


def load(name):
    spec = importlib.util.spec_from_file_location(name, REPO / f"code/requirements/{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


downloader = load("download_data")
installer = load("install_environments")
from file_io import directory_hashes
from runtime import sha256


class ResourcePreparationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.work = self.root / "downloads"
        self.work.mkdir()
        self.producer = self.root / "code/gse274113/prepare_scenicplus_offline_genome_resources.py"
        self.producer.parent.mkdir(parents=True)
        shutil.copy2(REPO / self.producer.relative_to(self.root), self.producer)
        self.producer_hash = sha256(self.producer)
        with gzip.open(self.work / "Homo_sapiens.GRCh38.116.gtf.gz", "wt") as handle:
            for index in range(15000):
                attrs = f'gene_id "G{index}"; gene_name "G{index}"; gene_biotype "protein_coding";'
                for feature in ("gene", "transcript"):
                    suffix = ' transcript_biotype "protein_coding";' if feature == "transcript" else ""
                    handle.write(f'1\tfixture\t{feature}\t1\t10\t.\t+\t.\t{attrs}{suffix}\n')
        self.sizes = self.work / "hg38.p14.chrom.sizes"
        self.chromosomes = "".join(f"chr{name}\t100\n" for name in [*range(1, 23), "X", "Y", "M"])
        self.sizes.write_text(self.chromosomes)
        self.output = self.root / "data/gse274113/external_resources/SCENICplus_GRCh38p14_Ensembl116"
        settings = patch.multiple(downloader, REPO=self.root, WORK=self.work, USER_HOME=self.root / "home",
                                  python_executable=lambda: sys.executable)
        settings.start()
        self.addCleanup(settings.stop)

    def prepare(self, force=False):
        task = downloader.Downloader(force=force)
        with patch.object(task, "fetch"), contextlib.redirect_stdout(io.StringIO()):
            task.scenicplus()

    def test_preparation_force_and_frozen_producer(self):
        self.prepare()
        first = directory_hashes(self.output)
        self.assertIn("artifact_manifest.tsv", first)
        self.prepare(force=True)
        self.assertEqual(directory_hashes(self.output), first)
        self.assertEqual(sha256(self.producer), self.producer_hash)
        backup = next((self.root / "home/backups/multiCSN_manuscript/scenicplus_resources").iterdir())
        self.assertEqual(directory_hashes(backup / "contents"), first)

    def test_failed_preparation_can_retry_and_preserves_existing_output(self):
        self.sizes.write_text("chr1\t100\n")
        with self.assertRaises(subprocess.CalledProcessError):
            self.prepare()
        self.assertFalse(self.output.exists())
        self.sizes.write_text(self.chromosomes)
        self.prepare()
        first = directory_hashes(self.output)
        self.sizes.write_text("chr1\t100\n")
        with self.assertRaises(subprocess.CalledProcessError):
            self.prepare(force=True)
        self.assertEqual(directory_hashes(self.output), first)
        self.assertFalse(list(self.output.parent.glob(".scenicplus-resources-*")))


class LingerInstallationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.prefix = self.root / "envs/linger"
        venv.EnvBuilder(with_pip=False, symlinks=True).create(self.prefix)
        python = self.prefix / "bin/python"
        site = Path(subprocess.check_output([str(python), "-c", "import sysconfig;print(sysconfig.get_path('purelib'))"], text=True).strip())
        self.package = site / "LingerGRN"
        self.package.mkdir()
        metadata = site / "LingerGRN-1.110.dist-info"
        metadata.mkdir()
        (metadata / "METADATA").write_text("Metadata-Version: 2.1\nName: LingerGRN\nVersion: 1.110\n")
        (self.root / "sources/LINGER/.git").mkdir(parents=True)
        self.ll_net = self.package / "LL_net.py"
        self.ll_net.write_text('''def fixture():
    for ii in []:
        if True:
            TFName_temp=TFName[np.array(TFidxtemp)]
            temps=np.abs(AA0).mean(axis=0)
            zscored_arr = np.nan_to_num(temps, nan=0.0)
        for k in range(len(TFidxtemp)):
            TG_1.append(geneName[ii])
            TF_1.append(TFName_temp[k])
            score_1.append(zscored_arr[k])
    TF_TG=pd.DataFrame(TG_1)
    TF_TG.columns=['TG']
    TF_TG['TF']=TF_1
''')
        for name in ("LINGER_tr", "preprocess", "pseudo_bulk"):
            (self.package / f"{name}.py").write_text("")
        settings = patch.multiple(installer, ENV_ROOT=self.root / "envs", SOURCE_ROOT=self.root / "sources")
        settings.start()
        self.addCleanup(settings.stop)

    def test_namespace_package_is_patched_and_repeat_is_idempotent(self):
        task = installer.Installer()
        with patch.object(task, "run", wraps=task.run) as commands:
            task.linger()
            self.assertTrue(any("LingerGRN.LL_net" in str(call) for call in commands.call_args_list))
            self.assertFalse(any("pip" in str(call) for call in commands.call_args_list))
        self.assertIn("            for k in range(len(TFidxtemp)):", self.ll_net.read_text())
        first = self.ll_net.read_bytes()
        task.linger()
        self.assertEqual(self.ll_net.read_bytes(), first)

    def test_existing_environment_with_broken_submodule_is_rejected(self):
        (self.package / "pseudo_bulk.py").write_text("raise RuntimeError('broken submodule')\n")
        with self.assertRaises(subprocess.CalledProcessError):
            installer.Installer().linger()


if __name__ == "__main__":
    unittest.main()
