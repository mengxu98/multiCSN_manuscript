#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source code/functions/runtime.sh
run_task code/scmultisim/run_multiseed.R
run_task code/scmultisim/analyze_multiseed.R
run_task code/scmultisim/run_tf_gene_ablation.R
