"""
06_event_rp.py: build the event-return-period lookup table.

Reads:  config.yaml (glofas.*, simulation.n_years, output.rp_floor), files.sim_probs,
        files.events, files.sites
Writes: files.events_sim (event_id, sim_year, source, n_sites_rp_ge_2, n_sites_rp_ge_10),
        files.event_rp (event_id, site_id, p_event, aep, rp_years), only rows
        with rp_years >= output.rp_floor

Method:
- lambda = number of historical events / number of record years
  (glofas.end_year - glofas.start_year + 1).
- aep = 1 - exp(-lambda * (1 - p_event)); rp_years = 1 / aep.
- n_sites_rp_ge_2 and n_sites_rp_ge_10 count sites with rp_years >= 2 and >= 10
  over all of sim_probs, not the filtered event_rp. 2 and 10 are fixed column
  definitions, not config.
- The number of simulated years is simulation.n_years (years with no events are not
  in sim_probs); stop if any sim_year is above it.
- Print the number of simulated years, events and rows in event_rp; the largest
  rp_years overall, at each site (median, min and max across sites) and at the
  required sites; and the distribution of n_sites_rp_ge_10 per event.
"""
import numpy as np
import pandas as pd
import yaml

cfg = yaml.safe_load(open("config.yaml"))
f = cfg["files"]
n_record = cfg["glofas"]["end_year"] - cfg["glofas"]["start_year"] + 1
n_years = cfg["simulation"]["n_years"]

events = pd.read_parquet(f["events"])
lam = len(events) / n_record
sim = pd.read_parquet(f["sim_probs"])
if sim["sim_year"].max() > n_years:
    raise ValueError(f"sim_year up to {sim['sim_year'].max()} but simulation.n_years = {n_years}")
if not sim["p_event"].between(0, 1, inclusive="neither").all():
    raise ValueError("p_event outside (0, 1)")
sim["aep"] = 1 - np.exp(-lam * (1 - sim["p_event"]))
sim["rp_years"] = 1 / sim["aep"]

events_sim = (sim.assign(ge_2=sim["rp_years"] >= 2, ge_10=sim["rp_years"] >= 10)
              .groupby("event_id")
              .agg(sim_year=("sim_year", "first"), source=("source", "first"),
                   n_sites_rp_ge_2=("ge_2", "sum"), n_sites_rp_ge_10=("ge_10", "sum"))
              .reset_index())
event_rp = sim.loc[sim["rp_years"] >= cfg["output"]["rp_floor"],
                   ["event_id", "site_id", "p_event", "aep", "rp_years"]]
events_sim.to_parquet(f["events_sim"], index=False)
event_rp.to_parquet(f["event_rp"], index=False)

pd.set_option("display.width", 200)
print(f"lambda = {len(events)} events / {n_record} years = {lam:.3f}")
print(f"{n_years} simulated years ({events_sim['sim_year'].nunique()} with events), {len(events_sim)} events")
print(f"wrote {f['events_sim']} ({len(events_sim)} rows) and {f['event_rp']} "
      f"({len(event_rp)} rows with rp_years >= {cfg['output']['rp_floor']} of {len(sim)})")

cols = ["site_id", "event_id", "sim_year", "source", "p_event", "rp_years"]
print("\nLargest rp_years overall:")
print(sim.loc[[sim["rp_years"].idxmax()], cols].round({"rp_years": 0}).to_string(index=False))

sites = pd.read_parquet(f["sites"])[["site_id", "name"]]
top = sim.loc[sim.groupby("site_id")["rp_years"].idxmax(), cols].merge(sites, on="site_id")
top = top.sort_values("rp_years", ascending=False)
print(f"\nLargest rp_years at each of the {len(top)} sites:")
print(top.round({"rp_years": 0}).to_string(index=False))
print(f"across sites: median {top['rp_years'].median():.0f}, min {top['rp_years'].min():.0f}, "
      f"max {top['rp_years'].max():.0f}")
print("\nRequired sites:")
print(top[top["name"] != ""].round({"rp_years": 0}).to_string(index=False))

n10 = events_sim["n_sites_rp_ge_10"]
print(f"\nn_sites_rp_ge_10 per event: median {n10.median():.0f}, 95th {n10.quantile(0.95):.1f}, "
      f"99th {n10.quantile(0.99):.1f}, max {n10.max()}")
print(f"n_sites_rp_ge_2 per event: median {events_sim['n_sites_rp_ge_2'].median():.0f}, "
      f"max {events_sim['n_sites_rp_ge_2'].max()}")
