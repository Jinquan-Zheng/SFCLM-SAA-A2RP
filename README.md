# SFCLM: Stochastic Flow-Capturing Location Model

This repository provides open-source code for solving a **two-stage stochastic flow-capturing location model (SFCLM)** for strategic transfer location planning in multi-hop ride-pooling systems. The SAA problem is solved with **Gurobi**. Solution quality is assessed by the **average two-replication procedure (A2RP)**, which supports **adaptive SAA sample size assessment** via the relative optimality gap (ROG).

> Accompanying paper:  A flow-capturing location model for multi-hop ride-pooling systems with stochastic demand. *Transportmetrica A: Transport Science*, 1–35. https://doi.org/10.1080/23249935.2026.2753573

## Citation

If you use this code, please cite:

Zheng, J., Li, B., Szeto, W. Y., Lei, T., Gong, L., & Luo, Q. (2026). A flow-capturing location model for multi-hop ride-pooling systems with stochastic demand. *Transportmetrica A: Transport Science*, 1–35. https://doi.org/10.1080/23249935.2026.2753573

```bibtex
@article{zheng2026sfclm,
  author  = {Zheng, J. and Li, B. and Szeto, W. Y. and Lei, T. and Gong, L. and Luo, Q.},
  title   = {A flow-capturing location model for multi-hop ride-pooling systems with stochastic demand},
  journal = {Transportmetrica A: Transport Science},
  year    = {2026},
  pages   = {1--35},
  doi     = {10.1080/23249935.2026.2753573}
}
```

## What is included

- Two-stage SFCLM formulation solved by Gurobi
- A2RP evaluation and ROG-based sample size adequacy check
- Bundled a subset of **sample scenarios** and **2019 Shenzhen metro candidate stations** for a runnable demo

## Requirements

- Python 3.11+ (recommended)
- A **Jupyter** environment (`jupyter` / `notebook` / `ipykernel`) to open and run `notebooks/SFCLM_gurobi_a2rp.ipynb`
- A valid **Gurobi** license (academic or commercial), with `gurobipy` available in the **same Python kernel** used by the notebook
- Python packages listed in `requirements.txt` (`pandas`, `numpy`, `scipy`, `tqdm`, `gurobipy`, and Jupyter-related packages)

Install dependencies in the conda/virtual environment that you will use as the notebook kernel:

```bash
pip install -r requirements.txt
```

## Quick start

1. Launch Jupyter (`jupyter notebook` or `jupyter lab`) and open `notebooks/SFCLM_gurobi_a2rp.ipynb`.
2. Select a kernel that has Gurobi and the packages above installed.
3. Run **Cell 1** (core functions), then **Cell 2** (demo main).
4. Results are written to `results/`.

Demo defaults (`DEMO_MODE = True` in Cell 2):

- Scenario folder: `data/sample_scenarios/` (15 scenarios)
- Facility counts: `p = 1..20`
- ROG tolerance: `xi = 0.001`

Alternatively, from the repository root (without opening Jupyter):

```bash
python run_demo.py
```

## Paper-scale experiments

To reproduce the paper's sample-size assessment:

1. Prepare SAA scenario folders for increasing sample sizes `M`, each containing scenario CSV files with the schema in `data/README_data.md`.
2. Point Cell 2 to those folders and set `DEMO_MODE = False` (then `p = 1..191`).
3. For each sample size `M`, solve the SFCLM and compute ROG by A2RP; increase `M` until ROG is below `xi` for all `p` values.

Full scenario data are available from the corresponding author upon reasonable request.

## Repository layout

```
sfclm-saa/
├── notebooks/SFCLM_gurobi_a2rp.ipynb
├── data/
│   ├── candidate_stations_2019.csv
│   ├── sample_scenarios/
│   └── README_data.md
├── results/
├── run_demo.py
├── requirements.txt
└── LICENSE
```


## License

MIT License. See `LICENSE`.

## Disclaimer

Sample trip scenarios are provided for **method demonstration only**. They are derived from historical taxi OD records after mode-choice simulation and path generation. They are not an official open traffic dataset release.
