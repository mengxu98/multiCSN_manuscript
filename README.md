# multiCSN manuscript

Analysis, figures and Source Data for the multiCSN manuscript.

Study archive: [Figshare](https://doi.org/10.6084/m9.figshare.34754148).

Run commands from the repository root.

To run the full analyses, install dependencies and prepare the additional raw data:

Activate the analysis Python and R environments first. General tasks use the current Python and `PATH` R. Existing BEELINE Conda environments (`beeline_py39`, `beeline_r453`), the GRNBoost2 environment at `~/.venvs/beeline-unified-py39`, and Pando's `../comparator_r_library/` are reused when present. LINGER, SCENIC+ and new GRNBoost2 installations use `../comparator_envs/`. Conda is located through `CONDA_EXE` or `PATH`; GO resource preparation requires Bioconductor 3.22.

```bash
bash code/requirements/install_environments.sh
bash code/requirements/download_data.sh
```

The installer updates multiCSN from GitHub while reusing the other installed packages. To use changes from a local multiCSN checkout, install it afterwards in the same R environment:

```bash
R CMD INSTALL ../multiCSN
```

Run the required tasks:

| Command | Task |
| --- | --- |
| `bash 01_beeline.sh` | BEELINE benchmark and scaling |
| `bash 02_scmultisim.sh` | Simulations and ablation |
| `bash 03_gse274113.sh` | Perturbation benchmark |
| `bash 04_gse162170.sh` | Cortex analysis |
| `bash 05_gse194122.sh` | Erythroid differentiation analysis |
| `bash 06_plotting.sh` | Figures |
| `Rscript --vanilla code/plotting/assemble_perturbation_validation.R` | Fig. 3 and Supplementary Figs. S8 and S9 |
| `Rscript --vanilla code/plotting/figs2_compact.R` | Supplementary Fig. S2 |
| `bash 07_source_data.sh` | Source Data |
| `python3 code/gse274113/run_timing_benchmark.py` | Optional matched-host timing |

To compare static and supplied-pseudotime inference on the existing 500 BEELINE cases, use a new output directory:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 LC_ALL=C Rscript --vanilla code/beeline/compare_static_pseudotime.R runs/beeline/static_comparison_new
```

This comparison requires the existing AP-inclusive case table and historical multiCSN networks under `results/beeline/`. It verifies all 500 historical network replays, runs five alternating inference-and-export timing repetitions per mode, and writes separate mode metrics and a joint comparison with the ten existing comparators. It preserves the formal benchmark. Current-host mode timings are reported separately from historical comparator timings.

The revised BEELINE figures use a corrected LEAP rerun and an isolated plotting view. The frozen rerun, shuffle networks, and full static/pseudotime runs are retained on the study compute server; `results/beeline/` contains the historical benchmark. The comparison and state-order control tables are under `source_data/supplementary_tables/`.

The comparison uses [compare_static_pseudotime.R](code/beeline/compare_static_pseudotime.R), [rerun_leap_ordered.R](code/beeline/rerun_leap_ordered.R), [run_pseudotime_shuffle.R](code/beeline/run_pseudotime_shuffle.R) and [assemble_pseudotime_revision.py](code/beeline/assemble_pseudotime_revision.py). Run `Rscript --vanilla code/beeline/test_leap_pseudotime.R` to verify the LEAP ordering. The corrected LEAP runtime is reported separately from archived matched timings. The 1,500 state-order controls are in [state_order_controls.tsv](source_data/supplementary_tables/beeline_pseudotime_20261105/state_order_controls.tsv).

Full analyses require Linux and the raw inputs declared by the analysis scripts.
Task parameters are defined in those scripts; preparation tasks build the required candidates and resources.

Original analysis and plotting code is licensed under the [MIT License](LICENSE).
Author-generated result tables, Source Data and figures are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Bundled third-party
code and reused data retain their original licenses; existing notices and
licenses under `code/vendor/` continue to apply.
