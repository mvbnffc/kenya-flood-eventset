# Kenya flood event set

Stochastic river-flood event set for Kenya, built from GloFAS v5 daily discharge with the Heffernan–Tawn conditional extremes model. Produces an event–return-period lookup table for the UNDRR Kenya disaster risk profile.

## Setup

1. `conda env create -f environment.yml`, then `conda activate kenya-flood-eventset`
2. In R: `install.packages(c("texmex", "arrow", "yaml"))`
3. Put your EWDS API key in `~/.cdsapirc` (see the "How to use the API" page on the Copernicus Early Warning Data Store).
4. Download these into `data/raw/` by hand from the GloFAS v5 auxiliary data page: the upstream-area map, the flow-direction (LDD) map and the return-period threshold file. File names are set in `config.yaml`.

## Run

From the repo root, in order:

```
python scripts/01_download.py
python scripts/02_sites.py
python scripts/03_extract.py
python scripts/04_events.py
Rscript scripts/05_simulate.R
python scripts/06_event_rp.py
python scripts/07_checks.py
```

All parameters are in `config.yaml`. The default config runs the small end-to-end version (80 sites, 1,000 years).

## Data

`data/`, `results/` and `figures/` are not tracked. GloFAS data come from the Copernicus Emergency Management Service; cite them as its licence requires.
