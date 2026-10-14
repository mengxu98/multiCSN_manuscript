# multiCSN manuscript

Analysis, figures and Source Data for the multiCSN manuscript.

Run commands from the repository root. Install dependencies and prepare data:

```bash
bash code/requirements/install_environments.sh
bash code/requirements/download_data.sh
```

The default package installer follows the development branches. To reproduce the
multiomics performance snapshot, extract `source_data/performance_reproducibility.zip`
and install its hash-verified source archives into a separate library:

```bash
Rscript code/requirements/install_csn_packages.R /path/to/extracted_snapshot /path/to/isolated_library
```

The snapshot manifest identifies the exact source files even when package versions
are unchanged. Add that library to the R library search path for the benchmark.

Run the required tasks:

| Command | Task |
| --- | --- |
| `bash 01_beeline.sh` | BEELINE benchmark and scaling |
| `bash 02_scmultisim.sh` | Simulations and ablation |
| `bash 03_gse274113.sh` | Perturbation benchmark |
| `bash 04_gse162170.sh` | Cortex analysis |
| `bash 05_gse194122.sh` | Erythroid differentiation analysis |
| `bash 06_plotting.sh` | Figures |
| `bash 07_source_data.sh` | Source Data |

Full analyses require Linux. Task parameters are defined in the analysis scripts.
