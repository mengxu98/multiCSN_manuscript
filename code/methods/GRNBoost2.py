"""Official Arboreto GRNBoost2 inference and network export."""
from pathlib import Path


def infer_tf_names(expression, path):
    """Preserve the candidate rule used in the existing BEELINE wrapper."""
    import pandas as pd
    genes = list(expression.columns)
    prefixed = [gene for gene in genes if gene.startswith('TF')]
    if prefixed:
        return prefixed
    metadata_path = path.with_name('case_metadata.csv')
    if metadata_path.exists():
        metadata = pd.read_csv(metadata_path)
        if not metadata.empty and 'RequestedTFs' in metadata:
            count = int(metadata.iloc[0]['RequestedTFs'])
            if count > 0:
                return genes[:min(count, len(genes))]
    return genes


def infer_grnboost2(expression, regulators, out_file, client, seed=1, verbose=False):
    from arboreto.algo import grnboost2

    network = grnboost2(expression_data=expression, tf_names=regulators,
                       client_or_address=client, seed=seed, verbose=verbose)
    Path(out_file).parent.mkdir(parents=True, exist_ok=True)
    table = network.rename(columns={"TF": "regulator", "importance": "weight"})[
        ["regulator", "target", "weight"]]
    table.to_csv(out_file, index=False)
    return table
