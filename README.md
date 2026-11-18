# Kenya flood event set

Stochastic river-flood event set for Kenya, built from GloFAS v5 daily discharge with the Heffernan–Tawn conditional extremes model. Produces an event–return-period lookup table for the UNDRR Kenya disaster risk profile.

## Setup

1. `conda env create -f environment.yml`, then `conda activate kenya-flood-eventset`
2. In R: `install.packages(c("texmex", "arrow", "yaml"))`. Needs texmex 2.4.9 (uses internal mexTransform).
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

## Method and references

Conditional extremes model of Heffernan and Tawn (2004, JRSS B), fitted with texmex. The event-set approach follows Keef, Tawn and Lamb (2013, Environmetrics) and Quinn et al. (2019, WRR, doi:10.1029/2018WR024205). Using modelled discharge rather than gauges follows Wing et al. (2020, WRR, doi:10.1029/2020WR027692) and Olcese et al. (2022, WRR, doi:10.1029/2022WR032743).

Differences from Quinn et al. (2019): dependence is fitted on event peaks (maximum in each event window), not on lagged 3-day-maximum series in a fixed 9-day window; there are no spatial limits (one national model, which fits inside the 1,500 km radius Olcese et al. used for Southern Africa).

Static maps: LISFLOOD static and parameter maps for GloFAS (Choulga et al. 2023, JRC, CC BY 4.0). Kenya boundary: geoBoundaries ADM0 (CC BY 4.0).

## Defaults to revisit after the end-to-end run

- Site density: `min_upstream_area_km2` and `area_growth_for_new_site`.
- Dry sites: `qc.max_zero_share`; dry-site floodplain zones need a return period from a nearby modelled site or a rainfall-based method.
- Event definition: trigger rate, merge gap, padding, and windows with no valley that stay longer than `split_above_days`. Test a fixed 9-day window (Quinn et al. 2019) as a sensitivity run through 04-07.
- Margin and dependence thresholds (0.8 / 0.8): check threshold-stability plots in 05.
- GPD shape: ξ < 0 at 22/50 sites gives hard upper limits (e.g. Ahero and site 76 only ~20% above the record maximum). Consider pooling ξ by region or constraining ξ ≥ 0; check the effect on the PML tail.
- Seasonal model: fit MAM and OND separately if dependence differs.
- Lower Athi/Galana: GloFAS dry-season zeros on a perennial river; check against WRA gauge data (Baricho).
- Check 5 in 07 waits for GloFAS v5 return-period thresholds.
- Fathom step: zones from HydroBASINS level 8 (level 10 where several sites share a unit), each taking its nearest site by distance and upstream area (Quinn et al. 2019).

