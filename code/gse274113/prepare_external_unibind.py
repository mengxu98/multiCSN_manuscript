from pathlib import Path
import csv, gzip, hashlib, json, zipfile

root = Path.cwd()
raw = root / 'data/gse274113/raw_external_downloads'
archive = raw / 'damo_hg38_TFBS_per_TF.zip'
with archive.open('rb') as handle:
    digest = hashlib.file_digest(handle, 'md5').hexdigest()
if digest != 'b6dd67a5bbaeb25c12b6ed45a935adf8':
    raise ValueError('UniBind download does not match the publisher MD5')
with (root / 'data/gse274113/heldout_pseudobulk_full/effect_model_summary.tsv').open() as handle:
    targets = {row['target'] for row in csv.DictReader(handle, delimiter='\t')}
assert len(targets) == 19
output = root / 'data/gse274113/external_resources/UniBind_perturbed_TFs_hg38.tsv.gz'
counts = []
with zipfile.ZipFile(archive) as z, gzip.open(output, 'wt') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['regulator', 'chrom', 'start', 'end'])
    for tf in sorted(targets):
        files = [n for n in z.namelist() if len(Path(n).parts) > 2 and Path(n).parts[1] == tf and n.endswith('.bed')]
        regions = set()
        input_rows = 0
        for name in files:
            for line in z.read(name).decode().splitlines():
                if not line or line.startswith(('#', 'track', 'browser')):
                    continue
                fields = line.split('\t')
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                assert start >= 0 and end > start, (name, line)
                regions.add((chrom, start + 1, end))
                input_rows += 1
        for region in sorted(regions):
            writer.writerow((tf, *region))
        counts.append(dict(regulator=tf, files=len(files), raw_rows=input_rows, unique_intervals=len(regions), source_files=files))
        print(tf, len(files), input_rows, len(regions), flush=True)
(raw / 'unibind_preparation.json').write_text(json.dumps({'record':'https://zenodo.org/records/15528243','version':'v1','release_date':'2025-05-27','md5':digest,'coordinates':'BED zero-based half-open converted to one-based closed; start + 1, end unchanged','selection':'All tissue contexts for each exact perturbed TF directory; duplicate intervals removed within TF','targets':counts}, indent=2))
