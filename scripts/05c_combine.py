"""
05c_combine.py: combine the cluster chunks into the simulated event pool.

Reads:  config.yaml (glofas.*, simulation.n_years, simulation.pool_factor, simulation.chunk_events),
        files.event_peaks, files.chunks/chunk_NNNN.parquet
Writes: files.sim_probs (event_id, source, site_id, p_event)

Method:
- n_pool and n_chunks as in 05_model.R: n_pool = round(pool_factor * lambda * n_years),
  lambda = rows of event_peaks / record years; n_chunks = ceil(n_pool / chunk_events).
- Stop if any of chunk_0001 ... chunk_<n_chunks> is missing or other chunk files are
  present, if any p_event is NA or outside (0, 1), if the event_ids are not exactly
  1 ... n_pool, or if any event lacks a site.
"""
import glob
import os

import numpy as np
import pandas as pd
import yaml

cfg = yaml.safe_load(open("config.yaml"))
f, sc = cfg["files"], cfg["simulation"]
n_record = cfg["glofas"]["end_year"] - cfg["glofas"]["start_year"] + 1
n_events = len(pd.read_parquet(f["event_peaks"]))
n_pool = round(sc["pool_factor"] * n_events / n_record * sc["n_years"])
n_chunks = int(np.ceil(n_pool / sc["chunk_events"]))

expected = {os.path.join(f["chunks"], f"chunk_{i:04d}.parquet") for i in range(1, n_chunks + 1)}
found = set(glob.glob(os.path.join(f["chunks"], "chunk_*.parquet")))
if expected - found:
    raise ValueError(f"{len(expected - found)} of {n_chunks} chunks missing, e.g. {sorted(expected - found)[:5]}")
if found - expected:
    raise ValueError(f"unexpected chunk files: {sorted(found - expected)[:5]}")
sim = pd.concat([pd.read_parquet(p) for p in sorted(expected)], ignore_index=True)
if sim["p_event"].isna().any():
    raise ValueError(f"{sim['p_event'].isna().sum()} NA p_event values")
if not sim["p_event"].between(0, 1, inclusive="neither").all():
    raise ValueError("p_event outside (0, 1)")
ids = np.sort(sim["event_id"].unique())
if not np.array_equal(ids, np.arange(1, n_pool + 1)):
    raise ValueError(f"event_ids are not 1 ... {n_pool}: {len(ids)} distinct, {ids.min()} to {ids.max()}")
n_sites = sim["site_id"].nunique()
if (sim.groupby("event_id").size() != n_sites).any():
    raise ValueError(f"not every event has all {n_sites} sites")
sim.to_parquet(f["sim_probs"], index=False)
print(f"wrote {f['sim_probs']}: {n_pool} events x {n_sites} sites from {n_chunks} chunks, "
      f"{(sim.drop_duplicates('event_id')['source'] == 'model').sum()} from the model")
