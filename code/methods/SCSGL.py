"""SCSGL inference and network export."""
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


def infer_scsgl(expression, out_file):
    import pandas as pd

    sys.path.insert(0, str(REPO / "code" / "vendor" / "SCSGL"))
    import run_scSGL

    Path(out_file).parent.mkdir(parents=True, exist_ok=True)
    scsgl_raw = Path(str(out_file) + ".scsgl_raw.tsv")
    run_scSGL.run_inference(expression, str(scsgl_raw))
    # Convert run_scSGL Gene1/Gene2/EdgeWeight output to regulator,target,weight CSV.
    raw = pd.read_csv(scsgl_raw, sep="\t")
    table = raw.rename(
        columns={"Gene1": "regulator", "Gene2": "target", "EdgeWeight": "weight"}
    )
    table = table[["regulator", "target", "weight"]].sort_values(
        "weight", ascending=False
    )
    table.to_csv(out_file, index=False)
    scsgl_raw.unlink()
