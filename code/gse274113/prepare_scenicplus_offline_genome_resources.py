#!/usr/bin/env python3
"""Build deterministic SCENIC+ GRCh38.p14 genome resources from local inputs."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
import re
import types


PRIMARY = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
UCSC = {name: ("chrM" if name == "MT" else f"chr{name}") for name in PRIMARY}
ANNOTATION_COLUMNS = [
    "Chromosome", "Start", "End", "Strand", "Gene",
    "Transcription_Start_Site", "Transcript_type",
]
CHROMSIZE_COLUMNS = ["Chromosome", "Start", "End"]


def parse_attributes(value: str) -> dict[str, str]:
    return dict(re.findall(r'([A-Za-z0-9_]+) "([^"]*)";', value))


def build_annotation(gtf_path: Path) -> list[dict[str, str | int]]:
    genes: dict[str, tuple[str, int, int, str, str]] = {}
    transcripts: list[tuple[str, str, int, int, str, str]] = []
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] not in UCSC:
                continue
            feature = fields[2]
            if feature not in {"gene", "transcript"}:
                continue
            attrs = parse_attributes(fields[8])
            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", "")
            biotype = attrs.get("gene_biotype", attrs.get("gene_type", ""))
            if not gene_id or not gene_name or biotype != "protein_coding":
                continue
            chrom = UCSC[fields[0]]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            if feature == "gene":
                genes[gene_id] = (chrom, start, end, strand, gene_name)
            else:
                transcript_biotype = attrs.get(
                    "transcript_biotype", attrs.get("transcript_type", "")
                )
                if transcript_biotype == "protein_coding":
                    transcripts.append(
                        (gene_id, chrom, start, end, strand, gene_name)
                    )
    rows: list[dict[str, str | int]] = []
    for gene_id, chrom, transcript_start, transcript_end, strand, gene_name in transcripts:
        if gene_id not in genes:
            continue
        gene_chrom, gene_start, gene_end, gene_strand, canonical_name = genes[gene_id]
        if chrom != gene_chrom or strand != gene_strand:
            raise RuntimeError(f"Inconsistent transcript/gene coordinates for {gene_id}")
        tss = transcript_start if strand == "+" else transcript_end
        rows.append({
            "Chromosome": chrom,
            "Start": gene_start,
            "End": gene_end,
            "Strand": strand,
            "Gene": canonical_name or gene_name,
            "Transcription_Start_Site": tss,
            "Transcript_type": "protein_coding",
        })
    unique = {tuple(row[column] for column in ANNOTATION_COLUMNS): row for row in rows}
    result = list(unique.values())
    result.sort(key=lambda row: (
        PRIMARY.index("MT" if row["Chromosome"] == "chrM" else str(row["Chromosome"])[3:]),
        int(row["Start"]), int(row["End"]), str(row["Gene"]),
        int(row["Transcription_Start_Site"]),
    ))
    return result


def build_chromsizes(path: Path) -> list[dict[str, int | str]]:
    wanted = set(UCSC.values())
    observed: dict[str, int] = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 2 and fields[0] in wanted:
                observed[fields[0]] = int(fields[1])
    missing = wanted - observed.keys()
    if missing:
        raise RuntimeError(f"Missing primary chromosomes: {sorted(missing)}")
    return [
        {"Chromosome": UCSC[name], "Start": 0, "End": observed[UCSC[name]]}
        for name in PRIMARY
    ]


def validate(
    annotation: list[dict[str, str | int]],
    chromsizes: list[dict[str, str | int]],
) -> None:
    if not annotation:
        raise RuntimeError("No protein-coding transcript annotations were produced")
    if len(chromsizes) != 25:
        raise RuntimeError(f"Expected 25 primary chromosomes, found {len(chromsizes)}")
    limits = {str(row["Chromosome"]): int(row["End"]) for row in chromsizes}
    genes = set()
    for row in annotation:
        chrom = str(row["Chromosome"])
        start, end = int(row["Start"]), int(row["End"])
        tss = int(row["Transcription_Start_Site"])
        if chrom not in limits or not (1 <= start <= tss <= end <= limits[chrom]):
            raise RuntimeError(f"Invalid annotation coordinate: {row}")
        if row["Strand"] not in {"+", "-"} or row["Transcript_type"] != "protein_coding":
            raise RuntimeError(f"Invalid annotation semantics: {row}")
        genes.add(str(row["Gene"]))
    if len(genes) < 15000:
        raise RuntimeError(f"Unexpectedly few protein-coding gene symbols: {len(genes)}")


def write_table(rows: list[dict[str, str | int]], columns: list[str], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# ---- Parameters ----
ARGUMENTS = {
    "gtf": Path("data/gse274113/external_resources/Homo_sapiens.GRCh38.116.gtf.gz"),
    "chrom_sizes": Path(
        "data/gse274113/external_resources/hg38.p14.chrom.sizes"
    ),
    "output": Path(
        "data/gse274113/external_resources/SCENICplus_GRCh38p14_Ensembl116"
    ),

}


def main() -> None:
    args = types.SimpleNamespace(**ARGUMENTS)
    args.output.mkdir(parents=True, exist_ok=False)

    annotation = build_annotation(args.gtf)
    chromsizes = build_chromsizes(args.chrom_sizes)
    validate(annotation, chromsizes)
    annotation_path = args.output / "genome_annotation.tsv"
    chromsizes_path = args.output / "chromsizes.tsv"
    write_table(annotation, ANNOTATION_COLUMNS, annotation_path)
    write_table(chromsizes, CHROMSIZE_COLUMNS, chromsizes_path)

    validation_rows = [
        {"field": "assembly", "value": "GRCh38.p14"},
        {"field": "annotation_release", "value": "Ensembl 116"},
        {"field": "chromosome_style", "value": "UCSC"},
        {"field": "primary_chromosomes", "value": str(len(chromsizes))},
        {"field": "annotation_rows", "value": str(len(annotation))},
        {"field": "unique_gene_symbols", "value": str(len({str(row['Gene']) for row in annotation}))},
    ]
    write_table(validation_rows, ["field", "value"], args.output / "resource_validation.tsv")
    print(f"annotation_rows\t{len(annotation)}")
    print(f"unique_genes\t{len({str(row['Gene']) for row in annotation})}")
    print(f"chromosomes\t{len(chromsizes)}")
    print("STATUS\tPASS")


if __name__ == "__main__":
    main()
