import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "functions"))
from runtime import REPO, use_python

if __name__ == "__main__":
    use_python("general")

import pandas as pd


def evaluate_locus_ordering(root=REPO, output=None):
    root = Path(root).resolve()
    output = Path(output).resolve() if output is not None else root / "results/gse162170/locus_display_ordering.tsv"
    screen=pd.read_csv(root/'results/gse162170/loci/all_target_locus_screen.tsv',sep='\t')
    current=screen[screen.eligible].sort_values(['sign_balance','median_chain_delta_bic','target_margin','chains','target'],ascending=[False,False,False,False,True]).copy()
    current['within_class_rank']=current.groupby('chromatin_anchor').cumcount()+1
    no_balance=current.sort_values(['median_chain_delta_bic','target_margin','chains','target'],ascending=[False,False,False,True]).copy()
    no_balance['no_balance_rank']=no_balance.groupby('chromatin_anchor').cumcount()+1
    no_sign=screen[screen.rna_matches_chromatin_anchor & screen.regions.between(5,15) & screen.chains.between(20,60) & (screen.median_chain_delta_bic>=3) & (screen.target_margin>=0.1)].copy()
    no_sign=no_sign.sort_values(['median_chain_delta_bic','target_margin','chains','target'],ascending=[False,False,False,True])
    no_sign['no_sign_filter_rank']=no_sign.groupby('chromatin_anchor').cumcount()+1
    r=current[current.target.isin(['ZNF804A','KCNC2'])][['target','chromatin_anchor','within_class_rank','positive','negative','median_chain_delta_bic']].merge(no_balance[['target','no_balance_rank']],on='target').merge(no_sign[['target','no_sign_filter_rank']],on='target')
    r.to_csv(output,sep='\t',index=False)


if __name__ == "__main__":
    evaluate_locus_ordering()
