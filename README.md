# Interconversion patterns and phenotype-specific inhibitors in neuroblastoma

Code and data for:

> **F. Covell** *[Mathematical model of phenotypic plasticity predicts interconversion is a drug target for neuroblastoma.]*. [2026].

A three-population ODE model of neuroblastoma – adrenergic (**N**), intermediate / I-type (**I**) and mesenchymal (**S**) cells – used to test how each of the 64 possible interconversion patterns changes the effect of an N- or S-targeting inhibitor.

This repository contains the files used to produce the paper.

---

## Layout

```
NB_Phenotypic_Plasticity/
├── 01_in_silico_drug_trial/     §2.3–2.4, §3.1–3.2 · Table 6, Figs 3–5, Tables A1–A2
│   ├── N_Drug/  relapse_sim_fixed_17-04-26.py      + Data/ (4 CSVs)
│   └── S_Drug/  MES_relapse_sim_fixed_17-04-26.py  + Data/ (4 CSVs)
├── 02_time_series_example/      Fig 2
│   ├── NdrugTS.py
│   └── SdrugTS.py
├── 03_stability/                §2.5, §3.3 · Figs 6–7
│   ├── Stability.ipynb
│   └── Data/ (3 CSVs)
├── 04_treatment_schedules/      §2.6, §3.4 · Table 4, Figs 8–11
│   ├── P1End23-06-2026.py            Schedule A – two-phase, no inhibitor
│   ├── SimpleSinhib18-06-2026.py     Schedule B – three-phase, S inhibitor
│   ├── SimpleNinhib18-06-2026.py     Schedule C – three-phase, N inhibitor
│   └── Complex4Cycle16-06-2026.py    Schedule D – four-phase, both inhibitors
└── analysis/
    ├── TSInsilico-24-08-26.R          Fig 2
    ├── Validation_Plots_18_06_2026.R  time-to-remission plots (run before FullAnalysis)
    └── FullAnalysis-28-07-2026.R      Table 6, Figs 3–11
```

## Which file makes which result

| Paper item | Simulation (Python) | Plotting (R) |
|---|---|---|
| Table 6, Tables A1–A2 | `01_in_silico_drug_trial/*` | `FullAnalysis` |
| Fig 2 | `02_time_series_example/NdrugTS.py`, `SdrugTS.py` | `TSInsilico` |
| Figs 3–5 (NTG and sensitivity) | `01_in_silico_drug_trial/*` | `FullAnalysis` |
| Figs 6–7 (stability) | `03_stability/Stability.ipynb` | `FullAnalysis` |
| Fig 8 (Schedule A) | `P1End23-06-2026.py` | `FullAnalysis` |
| Fig 9 (Schedule B) | `SimpleSinhib18-06-2026.py` | `FullAnalysis` |
| Fig 10 (Schedule C) | `SimpleNinhib18-06-2026.py` | `FullAnalysis` |
| Fig 11 (Schedule D) | `Complex4Cycle16-06-2026.py` | `FullAnalysis` |

The sensitivity analysis (Figs 4–5) does not use separate simulations: `FullAnalysis` splits the 36 runs per pattern from the drug trial by dose and by interconversion rate.

## The model

Each cell type has logistic growth, drug-induced death and interconversion (P = N + S + I, carrying capacity K = 1):

```
dN/dt = k1N·N(1−P) − k2N·N·D + k3SN·S + k3IN·I − k3NS·N − k3NI·N
dS/dt = k1S·S(1−P) − k2S·S·D + k3NS·N + k3IS·I − k3SN·S − k3SI·S
dI/dt = k1I·I(1−P) − k2I·I·D + k3SI·S + k3NI·N − k3IS·I − k3IN·I
```

| Parameter | Value |
|---|---|
| Initial cells (NTG) | 1 × 10⁻⁶ each |
| Growth rates k1 (N, I, S) | 0.75, 0.5, 0.25 |
| Target-cell death rate k2 | 1 |
| Non-target death rate k2 | {0, 0.5} |
| Drug dose D | {0.1, 0.5, 0.9} (set as `Ind / Deg` = {0.05, 0.25, 0.45} / 0.5) |
| Interconversion rate k3 | 0 for closed paths; {0.1, 0.5, 0.9} for open paths |

### Pattern notation

Each cell type sits between two brackets: the left bracket is its path to the previous cell in the cycle N → I → S → N, the right bracket its path to the next. `<` / `>` = open, `|` = closed. For example `<N> |I> |S|`: N converts to both S and I, I converts to S, S converts to nothing.

| ID | Pattern | Class | Paired with |
|---|---|---|---|
| P1 | `<N> \|I> \|S\|` | N-suppressing | S inhibitor |
| P2 | `<N> \|I><S\|` | N-suppressing | S inhibitor |
| P3 | `\|N\| <I> \|S>` | S-suppressing | N inhibitor |
| P4 | `\|N\| <I> \|S\|` | S-suppressing | N inhibitor |
| P5 | `\|N\| <I><S>` | S-suppressing | N inhibitor |

## How to run

**Requirements**
- Python 3 with `numpy`, `scipy`, `pandas`, `matplotlib`, `sympy` (for the stability notebook)
- R with `dplyr`, `tidyr`, `tidyverse`, `ggplot2`, `patchwork`, `scales`, `DBI`, `RSQLite`, `ggh4x`, `writexl`, `openxlsx`, `stringr`

**Paths:** the scripts still contain the original `C:/Users/...` paths. Edit the `save_dir` / `setwd()` / `db_path` lines at the top or bottom of each script to point at your own copy before running.

**Order**
1. `01_in_silico_drug_trial/N_Drug/relapse_sim_fixed_17-04-26.py` and `S_Drug/MES_relapse_sim_fixed_17-04-26.py`
2. `02_time_series_example/NdrugTS.py` and `SdrugTS.py`
3. `03_stability/Stability.ipynb`
4. The four scripts in `04_treatment_schedules/`
5. In one R session: `analysis/Validation_Plots_18_06_2026.R`, then `analysis/FullAnalysis-28-07-2026.R` (FullAnalysis reuses plots made by Validation_Plots)
6. `analysis/TSInsilico-24-08-26.R`

## Data in this folder

| File | Contents |
|---|---|
| `01_in_silico_drug_trial/*/Data/InitialParamSpace.csv` | parameter combinations |
| `…/BurnInFinal.csv` | populations at end of NTG, one row per run (64 patterns × 36) |
| `…/TreatmentFinal.csv` | populations at end of treatment |
| `…/RelapseFinal.csv` | populations at end of relapse |
| `03_stability/Data/stability_*FULL.csv` | eigenvalues at the origin across dose and k3 |

The per-run `.pkl` files and the SQLite time-series databases (`Time_Series.db`, `timeseries.db`) are regenerated by the scripts and are not stored here.

## Notes

- **`Complex4Cycle16-06-2026.py`** was reconstructed from the 3-cycle `Complex16-06-2026.py` by changing only the save folder (`3 cycle` → `4 cycle`) and `cycle_length` (3 → 4). The original 3-cycle version is not included here.
- **`FullAnalysis-28-07-2026.R`, line 812:** the two-phase plots (Fig 8) now read `4 cycle/P1End/timeseries.db`. An earlier version read the three-phase N-inhibitor database by mistake.

[![DOI](https://zenodo.org/badge/1437577015.svg)](https://doi.org/10.5281/zenodo.23744299)
## Licence

Released under the [MIT Licence](LICENSE).

## Contact

F. Covell – F.Covell@hull.ac.uk
Centre of Excellence for Data Science, AI and Modelling, University of Hull
