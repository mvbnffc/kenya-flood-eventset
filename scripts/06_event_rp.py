"""
06_event_rp.py: assemble simulated years from the event pool and build the event-return-period lookup table.

Reads:  config.yaml (glofas.*, simulation.n_years, simulation.count_model, simulation.class_rp,
        simulation.class_edges, simulation.region_persistence, simulation.regions,
        simulation.seed, output.rp_floor), files.sim_probs (the pool), files.hist_probs,
        files.events, files.sites
Writes: files.events_sim (event_id, sim_year, source, event_class, region, n_sites_rp_ge_2,
        n_sites_rp_ge_10), files.event_rp (event_id, sim_year, site_id, p_event, p_event_raw,
        aep, rp_years), only rows with rp_years >= output.rp_floor; delivered events only

Method:
- Simulated margins are remapped to ranks because Heffernan-Tawn simulation
  over-predicts exceedances at non-conditioning sites (r1 diagnostic: +27% per event
  at T = 2, +51% at T = 10). At each site, p = rank / (N + 1), ties averaged; 05's
  value is kept as p_event_raw. Stop unless every site has every pool event.
- Classes are assigned on pool ranks; delivered return periods are re-ranked over the
  assembled event set, so per-event margins hold by construction.
- Event class: the number of sites with rp_years >= simulation.class_rp, binned by the
  lower edges simulation.class_edges (class k = number of edges <= n; 0 small, 1 medium,
  2 large for [1, 5]). Site region: simulation.regions maps each site's outlet_id
  (files.sites) to a region; stop if a modelled site's outlet is not mapped. Event
  region: the region with the most sites at rp_years >= class_rp; for ties, and small
  events, the tied region holding the event's highest-RP site (rank ties broken by
  p_event_raw). Historical events use their hist_probs p_event
  and the same rp_years formula; pool events use the pool ranks.
- Assembly (simulation.count_model = empirical; stop otherwise): for each of the
  simulation.n_years years, resample one record year (zeros included) with replacement
  and take its events' (class, region) cells. Each such event is matched to a pool event
  from the same cell with probability simulation.region_persistence, otherwise from the
  same class in any region. Cell-matched draws go first, then class-only draws from what
  is left; both without replacement, in random order, with numpy's
  default_rng(simulation.seed). Stop if a cell or class runs out. Unused pool events are
  dropped. Before drawing, print historical, pool and strict-matching need per cell,
  with the pool events that strict matching (persistence 1) could never deliver
  (cells with no historical events); after, the delivered counts per cell. Spatial
  differences within a year are kept only as far as the four regions; season and
  climate state are not captured.
- Re-rank the delivered events at each site: p_event = rank / (N + 1), N = number of
  delivered events.
- Largest per-site rp_years is bounded by roughly the simulation length; this is
  intended, as the Fathom maps cap return periods well below that for n_years = 10,000.
- lambda = number of historical events / number of record years
  (glofas.end_year - glofas.start_year + 1).
- aep = 1 - exp(-lambda * (1 - p_event)), from the remapped p_event; rp_years = 1 / aep.
- n_sites_rp_ge_2 and n_sites_rp_ge_10 count sites with rp_years >= 2 and >= 10
  over all delivered events, not the filtered event_rp. 2 and 10 are fixed column
  definitions, not config. event_class and region are the pool-rank cell used in assembly.
- Print delivered counts by class against historical shares, how many delivered events
  change class after re-ranking, the number of simulated years, events and rows in
  event_rp, the largest rp_years overall and at the required sites, and the
  distribution of n_sites_rp_ge_10 per event.
"""
import numpy as np
import pandas as pd
import yaml
from scipy.stats import rankdata

cfg = yaml.safe_load(open("config.yaml"))
f, sc = cfg["files"], cfg["simulation"]
years = np.arange(cfg["glofas"]["start_year"], cfg["glofas"]["end_year"] + 1)
n_years = sc["n_years"]
if sc["count_model"] != "empirical":
    raise ValueError(f"simulation.count_model must be empirical, not {sc['count_model']}")
outlet = pd.read_parquet(f["sites"]).set_index("site_id")["outlet_id"]
site_region = outlet.map(sc["regions"])
regions = list(dict.fromkeys(sc["regions"].values()))


def to_rp(p):
    return 1 / (1 - np.exp(-lam * (1 - p)))


def rank_p(df):
    """Per-site rank of p_event_raw over the events in df, as p = rank / (N + 1)."""
    n = df["event_id"].nunique()
    return df.groupby("site_id")["p_event_raw"].transform(lambda p: rankdata(p) / (n + 1))


def event_class(df):
    n = df.assign(x=df["rp_years"] >= sc["class_rp"]).groupby("event_id")["x"].sum()
    return pd.Series(np.searchsorted(sc["class_edges"], n, side="right"), index=n.index)


def event_cell(df):
    """(class, region) per event; df has event_id, site_id, rp_years, p_event_raw."""
    d = df.assign(region=df["site_id"].map(site_region), x=df["rp_years"] >= sc["class_rp"])
    n = d.groupby(["event_id", "region"])["x"].sum().unstack()
    tied = n.eq(n.max(axis=1), axis=0).stack()
    d = d[tied.reindex(pd.MultiIndex.from_arrays([d["event_id"], d["region"]])).to_numpy()]
    top = d.sort_values(["rp_years", "p_event_raw"]).groupby("event_id").tail(1).set_index("event_id")["region"]
    return pd.DataFrame({"event_class": event_class(df), "region": top})


events = pd.read_parquet(f["events"])
lam = len(events) / len(years)
pool = pd.read_parquet(f["sim_probs"]).rename(columns={"p_event": "p_event_raw"})
if not pool["p_event_raw"].between(0, 1, inclusive="neither").all():
    raise ValueError("p_event outside (0, 1)")
n_pool = pool["event_id"].nunique()
if (pool.groupby("site_id").size() != n_pool).any():
    raise ValueError(f"not every site has all {n_pool} pool events")
missing = site_region.reindex(pool["site_id"].unique()).isna()
if missing.any():
    raise ValueError(f"simulation.regions has no region for the outlet of sites {list(missing.index[missing])}: "
                     f"outlet_ids {sorted(outlet[missing.index[missing]].unique())}")
pool["rp_years"] = to_rp(rank_p(pool))
pool_cell = event_cell(pool)

hist = pd.read_parquet(f["hist_probs"]).rename(columns={"p_event": "p_event_raw"})
hist["rp_years"] = to_rp(hist["p_event_raw"])
hist_cell = event_cell(hist)
classes = np.arange(len(sc["class_edges"]) + 1)
cells = pd.MultiIndex.from_product([classes, regions], names=["event_class", "region"])
year = events.set_index("event_id").loc[hist_cell.index, "start"].dt.year
year_counts = (pd.crosstab(year, [hist_cell["event_class"], hist_cell["region"]])
               .reindex(index=years, columns=cells, fill_value=0))

rng = np.random.default_rng(sc["seed"])
draw = year_counts.to_numpy()[rng.integers(0, len(years), n_years)]   # n_years x cells
cell_of = np.tile(np.arange(len(cells)), n_years).repeat(draw.ravel())  # one entry per needed event
sim_year = np.repeat(np.arange(1, n_years + 1), draw.sum(axis=1))
strict = rng.random(len(cell_of)) < sc["region_persistence"]

pd.set_option("display.width", 200)
tab = pd.DataFrame({"historical": hist_cell.value_counts(), "pool": pool_cell.value_counts()},
                   index=cells).fillna(0).astype(int)
tab["need, strict"] = draw.sum(axis=0)
tab["never deliverable, strict"] = np.where(tab["historical"] == 0, tab["pool"], 0)
tab[f"need, cell-matched at {sc['region_persistence']}"] = np.bincount(cell_of[strict], minlength=len(cells))
print(f"cells (event class x region), before drawing; {tab['never deliverable, strict'].sum()} pool events "
      "could never be delivered under strict matching:")
print(tab.to_string())

event_id = np.zeros(len(cell_of), dtype=int)
for k, (c, g) in enumerate(cells):
    ids = pool_cell.index[(pool_cell["event_class"] == c) & (pool_cell["region"] == g)].to_numpy()
    want = np.flatnonzero(strict & (cell_of == k))
    if len(want) > len(ids):
        raise ValueError(f"cell ({c}, {g}): {len(want)} events needed, {len(ids)} in the pool")
    event_id[want] = rng.permutation(ids)[:len(want)]
for c in classes:
    ids = np.setdiff1d(pool_cell.index[pool_cell["event_class"] == c], event_id)
    want = np.flatnonzero(~strict & (cells.get_level_values(0)[cell_of] == c))
    if len(want) > len(ids):
        raise ValueError(f"class {c}: {len(want)} class-only events needed, {len(ids)} left in the pool")
    event_id[want] = rng.permutation(ids)[:len(want)]
picked = pd.DataFrame({"event_id": event_id, "sim_year": sim_year}).join(pool_cell, on="event_id")
print("\ndelivered per cell:")
print(picked.groupby(["event_class", "region"]).size().reindex(cells, fill_value=0).unstack().to_string())

sim = pool.drop(columns="rp_years").merge(picked, on="event_id")
sim["p_event"] = rank_p(sim)
sim["aep"] = 1 - np.exp(-lam * (1 - sim["p_event"]))
sim["rp_years"] = 1 / sim["aep"]

events_sim = (sim.assign(ge_2=sim["rp_years"] >= 2, ge_10=sim["rp_years"] >= 10)
              .groupby("event_id")
              .agg(sim_year=("sim_year", "first"), source=("source", "first"), event_class=("event_class", "first"),
                   region=("region", "first"), n_sites_rp_ge_2=("ge_2", "sum"), n_sites_rp_ge_10=("ge_10", "sum"))
              .reset_index().sort_values(["sim_year", "event_id"]))
event_rp = sim.loc[sim["rp_years"] >= cfg["output"]["rp_floor"],
                   ["event_id", "sim_year", "site_id", "p_event", "p_event_raw", "aep", "rp_years"]]
events_sim.to_parquet(f["events_sim"], index=False)
event_rp.to_parquet(f["event_rp"], index=False)

print(f"\nlambda = {len(events)} events / {len(years)} years = {lam:.3f}")
final_class = event_class(sim).reindex(events_sim["event_id"]).to_numpy()
tab = pd.DataFrame({"historical": np.bincount(hist_cell["event_class"], minlength=len(classes)),
                    "pool": np.bincount(pool_cell["event_class"], minlength=len(classes)),
                    "delivered": np.bincount(events_sim["event_class"], minlength=len(classes)),
                    "delivered, re-ranked class": np.bincount(final_class, minlength=len(classes))}, index=classes)
tab.insert(1, "historical share", tab["historical"] / tab["historical"].sum())
tab.insert(4, "delivered share", tab["delivered"] / tab["delivered"].sum())
print(f"\nevent classes by sites at RP >= {sc['class_rp']}, lower edges {sc['class_edges']}:")
print(tab.round(3).to_string())
moved = pd.crosstab(events_sim["event_class"].to_numpy(), final_class, rownames=["pool class"],
                    colnames=["re-ranked class"])
print(f"{(events_sim['event_class'].to_numpy() != final_class).sum()} of {len(events_sim)} delivered events "
      "change class after re-ranking:")
print(moved.to_string())

print(f"\n{n_years} simulated years ({events_sim['sim_year'].nunique()} with events), {len(events_sim)} events "
      f"from a pool of {n_pool}")
print(f"wrote {f['events_sim']} ({len(events_sim)} rows) and {f['event_rp']} "
      f"({len(event_rp)} rows with rp_years >= {cfg['output']['rp_floor']} of {len(sim)})")

cols = ["site_id", "event_id", "sim_year", "source", "p_event", "rp_years"]
print("\nLargest rp_years overall:")
print(sim.loc[[sim["rp_years"].idxmax()], cols].round({"rp_years": 0}).to_string(index=False))
sites = pd.read_parquet(f["sites"])[["site_id", "name"]]
top = sim.loc[sim.groupby("site_id")["rp_years"].idxmax(), cols].merge(sites, on="site_id")
print("\nLargest rp_years at the required sites:")
print(top[top["name"] != ""].round({"rp_years": 0}).to_string(index=False))

n10 = events_sim["n_sites_rp_ge_10"]
print(f"\nn_sites_rp_ge_10 per event: median {n10.median():.0f}, 95th {n10.quantile(0.95):.1f}, "
      f"99th {n10.quantile(0.99):.1f}, max {n10.max()}")
print(f"n_sites_rp_ge_2 per event: median {events_sim['n_sites_rp_ge_2'].median():.0f}, "
      f"max {events_sim['n_sites_rp_ge_2'].max()}")
