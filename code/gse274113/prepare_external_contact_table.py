"""Extract both Table S2 loop sets; exclude the non-coordinate summary footer."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "functions"))
from runtime import use_python
if __name__ == "__main__":
    use_python("general")

import pandas as pd
root = Path('data/gse274113/raw_external_downloads')
parts = []
for sheet in ('Long_loops', 'HiCCUP_loops'):
    d = pd.read_excel(root / 'HSPC_Table_S2.xlsx', sheet_name=sheet)
    cols = ['anchor1_chrom','anchor1_start','anchor1_end','anchor2_chrom','anchor2_start','anchor2_end']
    d = d.iloc[:,:6].copy(); d.columns = cols
    d['source_row'] = d.index + 2
    d = d.dropna(subset=cols)
    for c in ['anchor1_start','anchor1_end','anchor2_start','anchor2_end']: d[c] = d[c].astype(int)
    for c in ['anchor1_chrom','anchor2_chrom']: d[c] = d[c].astype(str).map(lambda v: v if v.startswith('chr') else 'chr'+v)
    d['source_sheet'] = sheet
    d['loop_id'] = sheet + '_' + d['source_row'].astype(str)
    parts.append(d)
x = pd.concat(parts,ignore_index=True)
assert len(x) == 3091
x.to_csv(root / 'HSPC_loops_hg19.tsv',sep='\t',index=False)
