#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source code/functions/runtime.sh
run_task code/plotting/fig1.py
run_task code/plotting/fig2.R
run_task code/plotting/fig3.R
run_task code/plotting/figs7.R
run_task code/plotting/fig4.R
run_task code/plotting/figs8.R
run_task code/plotting/figs9.R
run_task code/plotting/fig5.R
run_task code/plotting/figs10.R
run_task code/plotting/figs1.R
run_task code/plotting/figs2.R
run_task code/plotting/figs3.R
run_task code/plotting/figs4.R
run_task code/plotting/figs5.R
run_task code/plotting/figs6.R
run_task code/plotting/figs11.R
run_task code/plotting/figs12.R
python3 code/plotting/verify_figures.py
