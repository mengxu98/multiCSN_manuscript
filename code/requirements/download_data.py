#!/usr/bin/env python3
"""Download or prepare declared inputs; verify checks the inventory."""
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import csv
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "functions"))
from runtime import REPO, USER_HOME, python_executable, rscript_executable, task_environment, sha256
from file_io import publish_directory

WORK = USER_HOME / "srcpkgs"
GTF_URL = "https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz"
CHROMSIZES_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p14/hg38.p14.chrom.sizes"
INVENTORY = [('data/beeline/Synthetic', 'BEELINE synthetic'), ('data/beeline/Curated', 'BEELINE curated'), ('data/gse162170/gse162170_candidates.rds', 'GSE162170 candidates'), ('data/gse194122/processed_multiome.h5ad', 'GSE194122 H5AD'), ('data/gse194122/erythroid_expression_subset.rds', 'fig5 expression subset'), ('data/gse274113/objects/GSE274113_control_full', 'GSE274113 control object'), ('data/gse274113/external_resources/primary_erythroid_ChIP_peaks_hg38.tsv.gz', 'external ChIP'), ('data/gse274113/external_resources/UniBind_perturbed_TFs_hg38.tsv.gz', 'external UniBind'), ('data/gse274113/external_resources/HSPC_loops_hg38.tsv.gz', 'external HSPC loops'), ('data/gse274113/external_resources/SCENICplus_GRCh38p14_Ensembl116/artifact_manifest.tsv', 'SCENIC+ offline resources'), ('data/gse274113/external_resources/SCENICplus_hg38_screen_v10_clust/artifact_manifest.tsv', 'SCENIC+ screen resources'), ('data/resources/LINGER/data_bulk/TFName.txt', 'LINGER reference data')]
EXTERNAL_URLS = {
    "GATA1_hg18.bed.gz": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1067nnn/GSM1067274/suppl/GSM1067274_Erythroid_GATA1_peaks.bed.gz",
    "NFE2_hg18.bed.gz": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1067nnn/GSM1067276/suppl/GSM1067276_Erythroid_NFE2_peaks.bed.gz",
    "hg18ToHg38.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz",
    "hg19ToHg38.over.chain.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    "damo_hg38_TFBS_per_TF.zip": "https://zenodo.org/records/15528243/files/damo_hg38_TFBS_per_TF.zip?download=1",
}
MOTIF_COMMIT = "8bcb5db8521e1acf0925452b2b49a1a961ba0bd0"
MOTIF_FILES = {
    "motifs.rda": "4de169756670c43116c21253dccf9d8cb358b4d84134809f454d9cf5c90f75a5",
    "motif2tf.rda": "b36331c17f754ce19e5759164bb3afb3bc332d23d9e458c77ffc309972cfcc3c",
    "motif2tf.tsv": "caa147d21ad6d745c520d8c5e2666b464d999f57b9a5993c3dd278ddd40146f5",
    "summary.tsv": "70af9d08db700a35ecc81946e8c10f76ab72437928362884c358ec377babea3b",
}
MOTIF_PROVENANCE = {
    "motif_collection": "Pando extended human motif collection",
    "motif_collection_commit": "3ffbc21bd827646476357acdc3a3a6ec9048052f",
    "motif_resource_commit": MOTIF_COMMIT,
    "motif_mapping_origin": "CIS-BP=711; JASPAR2020=1061; SEQ_SIMILARITY=575",
    "motif_count": "1590",
    "motif_mapping_tfs": "1372",
    "motif_mapping_links": "2347",
    "motif_p_cutoff": "5e-05",
    "motif_background": "subject",
    "motif_window": "7",
    "motif_score_mode": "binary",
    "motif_scan_settings_source": "Defaults in the historical motifmatchr call; checked in versions 1.28.0 and 1.32.0",
    "motif_collection_reference": "https://doi.org/10.1038/s41586-022-05279-8",
}
INVENTORY += [(f"data/resources/motifs/{name}", f"motif resource {name}")
              for name in (*MOTIF_FILES, "SHA256SUMS", "provenance.tsv")]


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
        for name, expected in MOTIF_FILES.items():
            self.fetch(f"https://raw.githubusercontent.com/mengxu98/datasets/{MOTIF_COMMIT}/motifs/{name}", destination / name)
            if not self.dry_run and sha256(destination / name) != expected:
                raise RuntimeError(f"Motif resource checksum mismatch: {name}")
        if not self.dry_run:
            (destination / "SHA256SUMS").write_text("".join(f"{digest}  {name}\n" for name, digest in MOTIF_FILES.items()))
            with (destination / "provenance.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(("field", "value"))
                writer.writerows(MOTIF_PROVENANCE.items())

    def beeline(self):
        destination = REPO / "data/beeline"
        if not (destination / "Synthetic").is_dir() or not (destination / "Curated").is_dir():
            source = WORK / "Beeline"
            self.run(["git", "clone", "-q", "--depth", "1", "https://github.com/Murali-group/Beeline.git", source])
            if self.dry_run:
                print(f"[dry-run] Copy {source / 'data'} into {destination}")
            else:
                shutil.copytree(source / "data", destination, dirs_exist_ok=True)
        if not self.dry_run:
            for group, expected in (("Synthetic", 300), ("Curated", 200)):
                count = len(list((destination / group).rglob("ExpressionData.csv")))
                print(f"{group}: {count} cases (expected {expected})")

    def gse274113(self):
        print("Raw GEO GSE274113 10x H5 and metadata must be in data/gse274113/GSE274113/.")
        print("Source: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274113/suppl/GSE274113_RAW.tar")
        self.run_r("build_control_object.R")
        self.run_r("build_heldout_truth.R")
        self.check("data/gse274113/objects/GSE274113_control_full", "control object")
        self.check("data/gse274113/heldout_pseudobulk_full", "heldout pseudobulk")

    def gse162170(self):
        if (REPO / "data/gse162170/gse162170_candidates.rds").is_file() and not self.force:
            print("GSE162170 candidates already exist")
        else:
            print("Place GEO GSE162170 raw matrices in data/gse162170/; construct paired candidates before 04_gse162170.sh.")

    def gse194122(self):
        if not (REPO / "data/gse194122/processed_multiome.h5ad").is_file():
            print("Place GEO GSE194122 processed_multiome.h5ad in data/gse194122/.")
        print("fig5.R generates the expression subset on first use.")
        self.check("data/gse194122/erythroid_expression_subset.rds", "figure 5 expression subset")

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
        if not (output / "artifact_manifest.tsv").is_file() or self.force:
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
                    f"chrom_sizes=Path({str(chromsizes)!r}), output=Path({str(candidate)!r}), "
                    f"gtf_url={GTF_URL!r}, chrom_sizes_url={CHROMSIZES_URL!r}); module['main']()")
                self.run([python_executable(), "-c", expression])
            if self.dry_run:
                prepare(output.with_name(output.name + ".staging"))
                print(f"[dry-run] Validate and publish resources to {output}; back up any previous directory")
            else:
                output.parent.mkdir(parents=True, exist_ok=True)
                with tempfile.TemporaryDirectory(prefix=".scenicplus-resources-", dir=output.parent) as temporary:
                    candidate = Path(temporary) / "resources"
                    prepare(candidate)
                    publish_directory(candidate, output,
                                      USER_HOME / "backups/multiCSN_manuscript/scenicplus_resources")
        print("Place screen rankings/scores feather files and motif annotations in data/gse274113/external_resources/SCENICplus_hg38_screen_v10_clust.")
        self.check("data/gse274113/external_resources/SCENICplus_hg38_screen_v10_clust/artifact_manifest.tsv", "screen resource manifest")

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
