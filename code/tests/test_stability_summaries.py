"""Rerun only stability reductions from saved networks, never inference."""
import os
import sys
from pathlib import Path
import subprocess
import tempfile
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'code/functions'))
from runtime import r_function_command
BASE = ROOT / 'results/gse274113/replicate_deletion_stability'
ENV = dict(os.environ, LC_ALL='C')


def main():
    with tempfile.TemporaryDirectory(prefix='stability-regression-') as directory:
        for script, group in [('summarize_leave_one_replicate.R', 'stability_summary'),
                              ('summarize_pairwise_fold_stability.R', 'pairwise_stability')]:
            output = Path(directory) / group
            command = r_function_command(
                'code/gse274113/' + script, script.removesuffix('.R'),
                timepoints=BASE / 'timepoints', output=output)
            subprocess.run(command, cwd=ROOT, env=ENV, check=True)
            expected = sorted(p.relative_to(BASE / group) for p in (BASE / group).rglob('*') if p.is_file())
            actual = sorted(p.relative_to(output) for p in output.rglob('*') if p.is_file())
            assert actual == expected, (group, actual, expected)
            for rel in expected:
                a = pd.read_csv(BASE / group / rel, sep='\t')
                b = pd.read_csv(output / rel, sep='\t')
                pd.testing.assert_frame_equal(a, b, check_dtype=False, rtol=1e-10, atol=1e-12)
            refused = subprocess.run(command, cwd=ROOT, env=ENV, capture_output=True, text=True)
            assert refused.returncode != 0 and 'Output already exists' in refused.stderr
            print(group, len(expected), 'tables match; overwrite refused', flush=True)


if __name__ == '__main__':
    main()
