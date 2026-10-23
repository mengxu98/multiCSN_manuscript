#!/usr/bin/env python3
"""Download or prepare declared inputs; verify checks the inventory."""
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "functions"))
from runtime import REPO, USER_HOME, python_executable, rscript_executable, task_environment
from file_io import publish_directory

WORK = USER_HOME / "srcpkgs"
GTF_URL = "https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz"
CHROMSIZES_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p14/hg38.p14.chrom.sizes"
INVENTORY = [('data/beeline/Synthetic', 'BEELINE synthetic'), ('data/beeline/Curated', 'BEELINE curated'), ('data/gse194122/processed_multiome.h5ad', 'GSE194122 H5AD'), ('data/gse274113/GSE274113', 'GSE274113 raw release'), ('data/gse274113/external_resources/primary_erythroid_ChIP_peaks_hg38.tsv.gz', 'external ChIP'), ('data/gse274113/external_resources/UniBind_perturbed_TFs_hg38.tsv.gz', 'external UniBind'), ('data/gse274113/external_resources/HSPC_loops_hg38.tsv.gz', 'external HSPC loops'), ('data/gse274113/external_resources/SCENICplus_GRCh38p14_Ensembl116/genome_annotation.tsv', 'SCENIC+ offline resources'), ('data/gse274113/external_resources/SCENICplus_hg38_screen_v10_clust/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather', 'SCENIC+ screen resources'), ('data/resources/LINGER/data_bulk/TFName.txt', 'LINGER reference data')]

EXTERNAL_URLS = {
    "GATA1_hg18.bed.gz": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1067nnn/GSM1067274/suppl/GSM1067274_Erythroid_GATA1_peaks.bed.gz",
    "NFE2_hg18.bed.gz": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1067nnn/GSM1067276/suppl/GSM1067276_Erythroid_NFE2_peaks.bed.gz",
    "hg18ToHg38.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz",
    "hg19ToHg38.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    "damo_hg38_TFBS_per_TF.zip": "https://zenodo.org/records/15528243/files/damo_hg38_TFBS_per_TF.zip?download=1",
}
MOTIF_COMMIT = "8bcb5db8521e1acf0925452b2b49a1a961ba0bd0"
MOTIF_FILES = ("motifs.rda", "motif2tf.rda")
PERTURBATION_REPLICATES = (*range(1, 11), 12, 13, 14, 16)
PERTURBATION_RELEASE_FILES = (
    "GSE274113_annotated_metadata.csv.gz",
    *(f"GSE274113_filtered_feature_bc_matrix_{replicate}.h5" for replicate in PERTURBATION_REPLICATES),
    *(f"GSE274113_rep{replicate}_filtered_feature_bc_matrix.h5" for replicate in PERTURBATION_REPLICATES),
    "supplementary/NIHMS2076308-supplement-Table_S2.xlsx",
    "supplementary/NIHMS2076308-supplement-Table_S3.xlsx",
)
INVENTORY += [(f"data/gse274113/GSE274113/{name}", f"GSE274113 release {name}") for name in PERTURBATION_RELEASE_FILES]
CORTEX_RELEASE_FILES = (
    "GSE162170_multiome_rna_counts.tsv.gz", "GSE162170_multiome_atac_counts.tsv.gz",
    "GSE162170_multiome_atac_consensus_peaks.txt.gz", "GSE162170_multiome_cell_metadata.txt.gz",
    "GSE162170_multiome_cluster_names.txt.gz",
)
INVENTORY += [(f"data/gse162170/release/{name}", f"GSE162170 paired release {name}") for name in CORTEX_RELEASE_FILES]
INVENTORY += [(f"data/resources/motifs/{name}", f"motif resource {name}")
              for name in MOTIF_FILES]


class Downloader:
    def __init__(self, dry_run=False, force=False):
        self.dry_run, self.force = dry_run, force

    def run(self, command, cwd=REPO):
        command = [str(value) for value in command]
        if self.dry_run:
            print("[dry-run] " + shlex.join(command))
        else:
            subprocess.run(command, cwd=cwd, env=task_environment(), check=True)

    def fetch(self, url, destination):
        if destination.is_file() and destination.stat().st_size and not self.force:
            print(f"Already exists: {destination}")
            return
        if not self.dry_run:
            destination.parent.mkdir(parents=True, exist_ok=True)
        self.run(["curl", "-fL", "--retry", "3", "--retry-delay", "5", "-o", destination, url])

    def check(self, relative, description):
        present = (REPO / relative).exists()
        print(f"{'Present' if present else 'Missing'}: {description} ({relative})")
        return present

    def run_r(self, script):
        self.run([rscript_executable(), "--vanilla", REPO / "code/gse274113" / script])

    def run_python(self, script):
        self.run([python_executable(), REPO / "code/gse274113" / script])

    def motifs(self):
        destination = REPO / "data/resources/motifs"
        for name in MOTIF_FILES:
            self.fetch(f"https://raw.githubusercontent.com/mengxu98/datasets/{MOTIF_COMMIT}/motifs/{name}", destination / name)

    def beeline(self):
        destination = REPO / "data/beeline"
        if not (destination / "Synthetic").is_dir() or not (destination / "Curated").is_dir():
            with tempfile.TemporaryDirectory(prefix="beeline-") as temporary:
                source = Path(temporary) / "Beeline"
                self.run(["git", "clone", "-q", "--depth", "1", "https://github.com/Murali-group/Beeline.git", source])
                if self.dry_run:
                    print(f"[dry-run] Copy {source / 'data'} into {destination}")
                else:
                    for group in ("Synthetic", "Curated"):
                        shutil.copytree(source / "data" / group, destination / group, dirs_exist_ok=True)
        if not self.dry_run:
            for group, expected in (("Synthetic", 300), ("Curated", 200)):
                count = len(list((destination / group).rglob("ExpressionData.csv")))
                print(f"{group}: {count} cases (expected {expected})")

    def gse274113(self):
        print("Raw GEO GSE274113 10x H5 and metadata must be in data/gse274113/GSE274113/.")
        print("Source: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274113/suppl/GSE274113_RAW.tar")
        self.check("data/gse274113/GSE274113", "raw release")

    def gse162170(self):
        for name in CORTEX_RELEASE_FILES:
            self.fetch(f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162170/suppl/{name}", REPO / "data/gse162170/release" / name)

    def gse194122(self):
        if not (REPO / "data/gse194122/processed_multiome.h5ad").is_file():
            print("Place GEO GSE194122 processed_multiome.h5ad in data/gse194122/.")


    def external(self):
        raw = REPO / "data/gse274113/raw_external_downloads"
        for name, url in EXTERNAL_URLS.items():
            self.fetch(url, raw / name)
        print("HSPC Table S2 requires publisher browser access; place HSPC_Table_S2.xlsx in", raw)
        for name in ("hg18ToHg38", "hg19ToHg38"):
            self.run(["gunzip", "-kf", raw / f"{name}.over.chain.gz"])
        self.run_r("prepare_external_chip.R")
        self.run_python("prepare_external_unibind.py")
        self.run_python("prepare_external_contact_table.py")
        self.run_r("prepare_external_contact.R")
        for relative, description in INVENTORY:
            if description in {"external ChIP", "external UniBind", "external HSPC loops"}:
                self.check(relative, description)

    def scenicplus(self):
        output = REPO / "data/gse274113/external_resources/SCENICplus_GRCh38p14_Ensembl116"
        if not (output / "genome_annotation.tsv").is_file() or self.force:
            gtf = WORK / "Homo_sapiens.GRCh38.116.gtf.gz"
            chromsizes = WORK / "hg38.p14.chrom.sizes"
            self.fetch(GTF_URL, gtf)
            self.fetch(CHROMSIZES_URL, chromsizes)
            producer = REPO / "code/gse274113/prepare_scenicplus_offline_genome_resources.py"
            def prepare(candidate):
                expression = (
                    "import runpy; from pathlib import Path; "
                    f"module=runpy.run_path({str(producer)!r}); "
                    f"module['ARGUMENTS'].update(gtf=Path({str(gtf)!r}), "
                    f"chrom_sizes=Path({str(chromsizes)!r}), output=Path({str(candidate)!r})); module['main']()")
                self.run([python_executable(), "-c", expression])
            if self.dry_run:
                prepare(output.with_name(output.name + ".staging"))
                print(f"[dry-run] Validate and replace resources at {output}")
            else:
                output.parent.mkdir(parents=True, exist_ok=True)
                with tempfile.TemporaryDirectory(prefix=".scenicplus-resources-", dir=output.parent) as temporary:
                    candidate = Path(temporary) / "resources"
                    prepare(candidate)
                    publish_directory(candidate, output)
        print("Place screen rankings/scores feather files and motif annotations in data/gse274113/external_resources/SCENICplus_hg38_screen_v10_clust.")
        self.check("data/gse274113/external_resources/SCENICplus_hg38_screen_v10_clust/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather", "screen motif rankings")

    def linger_bulk(self):
        print("Download and extract data_bulk.tar.gz from the Google Drive links in sources/LINGER/docs/{GRN_infer,Cancer}.md into data/resources/LINGER/.")
        self.check("data/resources/LINGER/data_bulk/TFName.txt", "LINGER reference data")

    def verify(self):
        present = [self.check(relative, description) for relative, description in INVENTORY]
        if not all(present):
            raise RuntimeError(f"Inventory incomplete: {sum(present)} present, {len(present) - sum(present)} missing")


SETTINGS = dict(parts=("all",), dry_run=False, force=False)


def main():
    selected = SETTINGS["parts"]
    names = ("motifs", "beeline", "gse274113", "gse162170", "gse194122", "external", "scenicplus", "linger-bulk", "verify")
    if set(selected) - set(names) - {"all"}:
        raise ValueError("Unknown dataset selection")
    downloader = Downloader(SETTINGS["dry_run"], SETTINGS["force"])
    for name in names:
        if "all" in selected or name in selected:
            print(f"Preparing {name}", flush=True)
            if name == "verify" and SETTINGS["dry_run"]:
                print("[dry-run] Verify the declared inventory after preparation")
            else:
                getattr(downloader, name.replace("-", "_"))()

if __name__ == "__main__":
    main()
