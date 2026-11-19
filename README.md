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

Conditional extremes model of Heffernan and Tawn (2004, JRSS B), fitted with texmex. The event-set approach follows Keef, Tawn and Lamb (2013, Environmetrics) and Quinn et al. (2019, WRR, doi:10.1029/2018WR024205). Using modelled discharge rather than gauges follows Wing et al. (2020, WRR, doi:10.1029/2020WR027692) and Olcese et al. (2022, WRR, doi:10.1029/2022WR032743). The planned pluvial extension and the Fathom loss step follow Olcese et al. (2024, WRR, doi:10.1029/2023WR036580).

Simulated margins are remapped to ranks in 06 (p = rank / (N + 1) per site over N events), because Heffernan–Tawn simulation over-predicts exceedances at non-conditioning sites (r1 diagnostic: +27% per event at T = 2, +51% at T = 10). The fitted dependence is kept; the largest per-site return period is bounded by roughly the simulation length, which is intended, as the Fathom maps cap return periods well below that for n_years = 10,000. 07 checks the delivered event_rp.

Simulated years are assembled by event class, with no extra fitting. 05 simulates a pool of `pool_factor` × λ × n_years events without years. 06 classes every event by its number of sites at RP ≥ 5 (small 0, medium 1–4, large ≥ 5; `class_rp`, `class_edges`) and by region (`regions`, four groups of sites by outlet: Lake Victoria, Turkana/Rift, Tana and upper Ewaso Ng'iro, Athi/coast), taking the region with the most sites at RP ≥ 5, or the region of the highest-RP site. Each simulated year then resamples one record year and, for each of its events, draws a pool event of the same class and region with probability `region_persistence` (0.8), otherwise of the same class in any region, without replacement; unused pool events are dropped. Classes are assigned on pool ranks; delivered return periods are re-ranked over the assembled event set, so per-event margins hold by construction. This keeps how many large events wet years have and, at the scale of the four regions, where they fall; within-year spatial differences by season or climate state (e.g. which basins a wet OND hits) are not captured otherwise.

Differences from Quinn et al. (2019): dependence is fitted on event peaks (maximum in each event window), not on lagged 3-day-maximum series in a fixed 9-day window; there are no spatial limits (one national model, which fits inside the 1,500 km radius Olcese et al. used for Southern Africa).

Static maps: LISFLOOD static and parameter maps for GloFAS (Choulga et al. 2023, JRC, CC BY 4.0). Kenya boundary: geoBoundaries ADM0 (CC BY 4.0).

## Defaults to revisit after the end-to-end run

- Site density: `min_upstream_area_km2` and `area_growth_for_new_site`.
- Dry sites: `qc.max_zero_share`; dry-site floodplain zones need a return period from a nearby modelled site or a rainfall-based method.
- Event definition: trigger rate, merge gap, padding, and windows with no valley that stay longer than `split_above_days`. Test a fixed 9-day window (Quinn et al. 2019) as a sensitivity run through 04-07.
- Margin and dependence thresholds (0.8 / 0.8): check threshold-stability plots in 05.
- GPD shape: ξ < 0 at 20/45 sites gives hard upper limits (e.g. Ahero only ~20% above the record maximum). Consider pooling ξ by region or constraining ξ ≥ 0; check the effect on the PML tail.
- Seasonal model: fit MAM and OND separately if dependence differs.
- Within-year repeat flooding at the same sites is under-represented (check 7: 0.42 vs 0.57 historical); finer regions or a different labelling rule may help if annual losses show it matters.
- Excluded: Ewaso Ng'iro (North) below the Lorian Swamp: transmission-loss reach where GloFAS flow does not grow with upstream area; unreliable for flood statistics; low exposure (`sites.excluded_sites`).
- Lower Athi/Galana: GloFAS dry-season zeros on a perennial river; check against WRA gauge data (Baricho).
- Check 5 in 07 waits for GloFAS v5 return-period thresholds.
- Event definition (Olcese et al. 2024): look for events with two separate hotspots in the viewer's random draws; they found these come from merging unrelated floods and fixed them with shorter windows. If present, try per-site event windows after the 9-day test.
- GPD thresholds: automatic per-site threshold choice (Murphy, Tawn & Varty 2023, as in Olcese et al. 2024) instead of a fixed 0.8, as a loop over texmex fits. Optional; may help sites such as Garissa where ξ is unstable above 0.8.
- Conditioning sites: the upper Ewaso Ng'iro (Archer's Post to Merti) and lower Athi lead the check 1 excess and conditioned the widest r1 event. If they still dominate, keep them as sites but don't let them condition events.
- Season labels: define hydrological seasons (long rains about March–July, short rains about October–January) before testing a seasonal split.
- OND underestimation: GloFAS shows OND 2023 as moderate on the Tana and almost dry on the Athi. Check against WRA gauges; look at Seven Forks reservoir behaviour in v5.
- Return periods are event-rate based (Poisson) while annual counts are overdispersed; this differs from annual-maximum AEP at T < 10. Confirm Fathom's return-period definition.
- Compute: `mexMonteCarlo` slows sharply with more sites. Time it at about 150 sites before the full network; plan parallel chunks (`parLapply`, `clusterSetRNGStream`) for 10,000 years.
- Climate-state regimes (possible extension): classify years by IOD/ENSO to capture which basins wet years hit.

## Fathom loss step (plan)

- Site placement for the full network: match sites to Fathom reach outlets, each linked to a GloFAS cell within 10 cells and 20% upstream area (Olcese et al. 2024). Ask Fathom whether its Kenya reach-outlet layer can be shared.
- Zones: each river site takes its HydroBASINS upstream area (level 8, level 10 where several sites share a unit). Rivers below the GloFAS threshold, down to Fathom's 50 km², take the return period of their nearest site in the same basin by distance and upstream area (Quinn et al. 2019).
- Small-river bracket: run losses twice, small-river units copying their site's return period (full dependence) and drawing their own (independent). AAL is the same; the PML gap shows how much the assumption matters.
- Pluvial extension: ERA5 daily rainfall at 0.25° tiles as extra sites in the same model (consistent with GloFAS forcing), each driving Fathom pluvial maps in its tile; merge fluvial and pluvial events that overlap in time and space. Priority: Nairobi, Mombasa, Kisumu. Check ERA5 extremes against KMD stations.
- No double counting: "clean" layers with pluvial flooding only where there is no fluvial flooding at that return period, and vice versa.
- Depths: exposure counted above 15 cm (fluvial) and 10 cm (pluvial), with 30 cm as a sensitivity test.
- Return periods between mapped levels are interpolated linearly; values above Fathom's largest return period take the top map.
- Validation: people affected against EM-DAT, sampling simulated periods of the record's length and keeping the largest events to match EM-DAT's coverage (Olcese et al. 2024). Use WorldPop for population.
- Message for partners: annual exposure at 1-in-100 is well below the population of the 100-year flood zone (about half in Olcese et al. 2024), because no single year floods every 100-year zone.

