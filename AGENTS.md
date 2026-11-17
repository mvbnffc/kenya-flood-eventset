# Kenya flood event set

Builds a stochastic river-flood event set for Kenya from GloFAS v5 daily discharge, using the Heffernan–Tawn conditional extremes model (R package texmex). The output is an event–return-period lookup table: for each simulated event, the return period at each river site. Losses are calculated elsewhere, not in this repo.

## Pipeline

| Script | Reads | Writes |
|---|---|---|
| 01_download.py | config | raw GloFAS files, one per year |
| 02_sites.py | upstream-area and flow-direction maps | sites table, site map |
| 03_extract.py | raw files, sites | daily discharge (date × site), site QC table |
| 04_events.py | discharge | events table, event peaks (event × site) |
| 05_simulate.R | event peaks, events | texmex fit, diagnostics, margins table, historical and simulated probabilities |
| 06_event_rp.py | simulated probabilities, events | simulated events table, event–RP lookup table |
| 07_checks.py | all outputs | check figures |

Exact paths and parameters are in config.yaml. Run every script from the repo root. Each script's header describes its method; follow it.

## Rules
- Plain numbered scripts, run in order. No packages, classes, CLI frameworks or plugin structures.
- All parameters and paths live in config.yaml. No hard-coded values, no other config layers.
- Allowed dependencies: xarray, netCDF4, pandas, pyarrow, numpy, scipy, matplotlib, pyyaml, cdsapi; R: texmex, arrow, yaml. Ask before adding anything.
- Use texmex for all extreme value fitting and simulation. Do not reimplement GPD fitting, Heffernan–Tawn or the simulation scheme in Python.
- Each script: short header saying inputs, outputs and purpose. Aim for under 200 lines.
- Write a function only if it is used more than once.
- Fail loudly. No broad try/except, no fallbacks, no silent defaults, no logging frameworks.
- Never modify data/raw/.
- No new files, tests, CI, Docker or type-checking unless asked.
- If a method choice is ambiguous, stop and ask. Do not build a more general solution "to be safe".
- When editing, prefer deleting code to adding it.
