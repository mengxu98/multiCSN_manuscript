# multiCSN manuscript

Analysis, figures and Source Data for the multiCSN manuscript.

Run commands from the repository root. Install dependencies and prepare data:

```bash
bash code/requirements/install_environments.sh
bash code/requirements/download_data.sh
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
| `bash 07_source_data.sh` | Source Data |
| `python3 code/gse274113/run_timing_benchmark.py` | Optional matched-host timing |

Full analyses require Linux and the raw inputs declared by the analysis scripts.
Task parameters are defined in those scripts; preparation tasks build the required candidates and resources.
