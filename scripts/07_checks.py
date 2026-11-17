"""
07_checks.py: check the simulated event set against the historical record.

Reads:  config.yaml, files.sites, files.hist_probs, files.sim_probs, files.events,
        files.events_sim, files.margins_params, files.rp_thresholds
Writes: figures/check_1_marginals.png ... figures/check_5_glofas_thresholds.png,
        and prints a one-line summary per check

Checks:
1. Marginals: at each site, the share of simulated years whose annual maximum RP
   is at least T should be about 1/T, for T = 2, 5, 10, 50.
2. Dependence: extremal dependence chi at p = 0.9 and 0.95 against distance between
   sites, historical vs simulated.
3. Event extent: number of sites with RP >= 2, 5 and 10 per event and per year,
   historical vs simulated.
4. Big years: where the 1997-98, 2019-20 and 2023-24 extents fall in the simulated
   distribution (they should be rare, not impossible).
5. GloFAS thresholds: discharge at RP 2, 5 and 20 implied by margins_params vs the
   GloFAS v5 threshold file. v5 thresholds include trend adjustments, so expect
   differences where site_qc shows a trend.
"""
import yaml

cfg = yaml.safe_load(open("config.yaml"))

# To be written in an agent session.
