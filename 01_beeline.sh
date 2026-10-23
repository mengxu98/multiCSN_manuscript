#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source code/functions/runtime.sh
python3 code/beeline/run_unified_benchmark.py
python3 code/beeline/promote_unified_batch.py
run_task code/beeline/ap_sensitivity.R
run_task code/beeline/prepare_scaling_data.R
run_task code/beeline/run_core_scaling.R
