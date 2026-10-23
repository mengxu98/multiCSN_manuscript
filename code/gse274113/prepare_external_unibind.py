from pathlib import Path
import csv, gzip, zipfile

root = Path.cwd()
raw = root / 'data/gse274113/raw_external_downloads'
archive = raw / 'damo_hg38_TFBS_per_TF.zip'
metadata = root / 'data/gse274113/GSE274113/GSE274113_annotated_metadata.csv.gz'
with gzip.open(metadata, 'rt') as handle:
    reader = csv.DictReader(handle)
    cell_columns = [name for name in reader.fieldnames if not name or name.startswith('Unnamed:')]
    if len(cell_columns) != 1 or 'target' not in reader.fieldnames:
        raise ValueError('Missing released cell identifiers or perturbation targets')
    targets = {row['target'] for row in reader
               if row[cell_columns[0]].strip() and row['target'] and row['target'] != 'NT'}
output = root / 'data/gse274113/external_resources/UniBind_perturbed_TFs_hg38.tsv.gz'
with zipfile.ZipFile(archive) as z, gzip.open(output, 'wt') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['regulator', 'chrom', 'start', 'end'])
    for tf in sorted(targets):
        files = [n for n in z.namelist() if len(Path(n).parts) > 2 and Path(n).parts[1] == tf and n.endswith('.bed')]
        regions = set()
        for name in files:
            for line in z.read(name).decode().splitlines():
                if not line or line.startswith(('#', 'track', 'browser')):
                    continue
                fields = line.split('\t')
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                regions.add((chrom, start + 1, end))
        for region in sorted(regions):
            writer.writerow((tf, *region))
        print(tf, len(files), len(regions), flush=True)
