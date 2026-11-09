# multiCSN manuscript

Analysis, figures and Source Data for the multiCSN manuscript.

Study archive: [Figshare](https://doi.org/10.6084/m9.figshare.34754148).

Run commands from the repository root.

To run the full analyses, activate the analysis Python and R environments, then install dependencies and prepare the raw data. The installer reuses existing method environments where available. Set `CONDA_EXE` if Conda is not on `PATH`. GO resource preparation requires Bioconductor 3.22.

```bash
bash code/requirements/install_environments.sh
bash code/requirements/download_data.sh
```

The installer updates multiCSN from GitHub while reusing the other installed packages. To use changes from a local multiCSN checkout, install it afterwards in the same R environment:

```bash
R CMD INSTALL ../multiCSN
```

Use the commands below for analysis, figures, and Source Data:

| Command | Task |
| --- | --- |
| `bash 01_beeline.sh` | BEELINE benchmark and scaling |
| `bash 02_scmultisim.sh` | Simulations and ablation |
| `bash 03_gse274113.sh` | Perturbation benchmark |
| `bash 04_gse162170.sh` | Cortex analysis |
| `bash 05_gse194122.sh` | Erythroid differentiation analysis |
| `bash 06_plotting.sh` | Run the main plotting workflow |
| `Rscript --vanilla code/plotting/figs2.R` | Supplementary Fig. S2 |
| `Rscript --vanilla code/plotting/figs8.R` | Fig. 3 and Supplementary Figs. S8 and S9 |
| `bash 07_source_data.sh` | Generate Source Data |
| `Rscript --vanilla code/beeline/compare_static_pseudotime.R runs/beeline/static_comparison_new` | Optional static and supplied-pseudotime comparison; requires prepared benchmark inputs and a new output directory. Set OpenMP and the active BLAS library to one thread for timing. |
| `python3 code/gse274113/run_timing_benchmark.py` | Optional matched-host timing |

Full analyses require Linux and the raw inputs declared by the analysis scripts.
Task parameters are defined in those scripts; preparation tasks build the required candidates and resources.

Original analysis and plotting code is licensed under the [MIT License](LICENSE).
Author-generated result tables, Source Data and figures are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Bundled third-party
code and reused data retain their original licenses; existing notices and
licenses under `code/vendor/` continue to apply.
