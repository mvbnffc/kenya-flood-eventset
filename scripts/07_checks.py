"""
07_checks.py: check the simulated event set against the historical record.

Reads:  config.yaml (glofas.*, simulation.n_years, output.rp_floor, checks.*), files.sites,
        files.hist_probs, files.event_rp, files.events, files.events_sim,
        files.sim_probs (check 6 only), files.dependence_params, files.boundary;
        files.margins_params and files.rp_thresholds for check 5
Writes: figures/check_1_marginals.png ... figures/check_6_event_profile.png
        (no check 5 figure while files.rp_thresholds is null; none for check 7), and prints
        a summary per check

This checks the delivered event_rp: simulated rp_years are read from files.event_rp
(after 06's assembly and re-ranking), with their sim_year. Event-site pairs missing
from event_rp are below output.rp_floor and count as rp_years 0, below every T; stop
if any T in checks.* is below output.rp_floor. Historical RP uses 06's formula:
rp_years = 1 / (1 - exp(-lambda * (1 - p_event))), lambda = historical events per record year.

Checks:
1. Marginals: at each site, the share of simulated years whose annual maximum RP
   is at least T, for T in checks.marginal_T, against the expected share given the
   historical annual counts N: 1 - mean over record years of (1 - q)^N, where
   q = -ln(1 - 1/T) / lambda is the per-event probability of RP >= T. Counts are
   overdispersed, so this is below 1/T at small T; 1/T is drawn as a reference line
   only. Print the per-site ratio of simulated to expected share for T in
   checks.ratio_T, sorted by the ratio at the second T, and its median, range and
   pooled value (all site-years). The expected share assumes independent events within
   a year, so the historical record's own ratio, printed alongside, is the reference.
2. Dependence: event co-occurrence rate, as in Wing et al. (2020) and Olcese et al.
   (2022). For each pair of sites i, j: COR(T) = share of events with RP >= T at i
   that also have RP >= T at j. Compute for T in checks.cor_T, historical vs
   simulated, plotted against great-circle distance between sites. Print the
   root-mean-square difference (CORE) over all pairs and over pairs with simulated
   COR > checks.cor_min, and the mean difference. Historical COR is noisy with 46
   years, so look for bias, not pair-by-pair agreement.
3. Event extent: number of sites with RP >= T (checks.extent_T) per event and per
   year, historical vs simulated, as exceedance curves.
4. Big years: where the extents of the calendar years in checks.big_years fall in
   the simulated annual distribution (they should be rare, not impossible).
5. GloFAS thresholds: discharge at RP 2, 5 and 20 implied by margins_params vs the
   GloFAS v5 threshold file. v5 thresholds include trend adjustments, so expect
   differences where site_qc shows a trend. Skipped while files.rp_thresholds is
   null; stops if it is set, as the check is not written yet.
6. Event profile: for simulated event checks.profile_event (must be a model event),
   the conditioning site (the site with the highest p_event in sim_probs, before
   06's remap, as mexMonteCarlo keeps only draws where the conditioning site is
   largest), a map of rp_years at all sites (hollow below output.rp_floor), and Heffernan-Tawn alpha and beta from the conditioning site to every
   other site against distance.
7. Repeat flooding: site-years with RP >= checks.repeat_T in two or more events, as a
   share of all site-years and of site-years with at least one, historical vs
   simulated; plus the median and range across sites of the second share.
"""
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.colors import LogNorm

cfg = yaml.safe_load(open("config.yaml"))
f, ck = cfg["files"], cfg["checks"]
years = np.arange(cfg["glofas"]["start_year"], cfg["glofas"]["end_year"] + 1)
sim_years = np.arange(1, cfg["simulation"]["n_years"] + 1)
HIST, SIM = "#404040", "#2a78b5"


def annual_max(df, index):
    """Annual maximum rp_years (year x site); years with no events, and site-years below rp_floor, get 0."""
    return df.groupby(["year", "site_id"])["rp_years"].max().unstack().reindex(index=index, fill_value=0)[site_ids].fillna(0)


def wide(df):
    return df.pivot(index="event_id", columns="site_id", values="rp_years")[site_ids]


def km(lon1, lat1, lon2, lat2):
    """Great-circle distance in km."""
    lon1, lat1, lon2, lat2 = map(np.radians, (lon1, lat1, lon2, lat2))
    h = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * 6371 * np.arcsin(np.sqrt(h))


def expected(T):
    """Expected share of years with annual max RP >= T, given the historical annual counts."""
    return 1 - np.mean((1 + np.log(1 - 1 / T) / lam) ** n_hist)


events = pd.read_parquet(f["events"])
lam = len(events) / len(years)
n_hist = events["start"].dt.year.value_counts().reindex(years, fill_value=0).values
hist = pd.read_parquet(f["hist_probs"]).merge(events[["event_id", "start"]], on="event_id")
hist["year"] = hist.pop("start").dt.year
hist["rp_years"] = 1 / (1 - np.exp(-lam * (1 - hist["p_event"])))
check_T = ck["marginal_T"] + ck["ratio_T"] + ck["cor_T"] + ck["extent_T"]
if min(check_T) < cfg["output"]["rp_floor"]:
    raise ValueError(f"check T = {min(check_T)} is below output.rp_floor = {cfg['output']['rp_floor']}")
events_sim = pd.read_parquet(f["events_sim"])
sim = pd.read_parquet(f["event_rp"]).rename(columns={"sim_year": "year"})
site_ids = np.sort(sim["site_id"].unique())
if not np.array_equal(np.sort(hist["site_id"].unique()), site_ids):
    raise ValueError("hist_probs and event_rp have different sites")
sites = pd.read_parquet(f["sites"]).set_index("site_id").loc[site_ids]
amax_h, amax_s = annual_max(hist, years), annual_max(sim, sim_years)
ev_h, ev_s = wide(hist), wide(sim).reindex(events_sim["event_id"]).fillna(0)
os.makedirs(f["figures"], exist_ok=True)
pd.set_option("display.width", 200)
print(f"lambda = {lam:.3f}; {len(site_ids)} sites; {len(ev_h)} historical events in {len(years)} years, "
      f"{len(ev_s)} simulated events in {len(sim_years)} years")

# Check 1: marginals.
fig, ax = plt.subplots(figsize=(7, 5))
Ts = np.array(ck["marginal_T"])
for T in Ts:
    share = (amax_s >= T).mean()
    ax.scatter(T * np.exp(np.random.default_rng(0).uniform(-0.08, 0.08, len(share))), share, s=10, color=SIM, alpha=0.6)
ax.plot(Ts, [expected(T) for T in Ts], "o-", color=HIST, label="expected from historical annual counts")
ax.plot(Ts, 1 / Ts, "--", color="0.6", label="1/T (reference)")
ax.scatter([], [], s=10, color=SIM, label="simulated, one point per site")
ax.set(xscale="log", yscale="log", xlabel="T (years)", ylabel="share of simulated years with annual max RP >= T",
       title="Check 1: marginals")
ax.set_xticks(Ts, [str(T) for T in Ts], minor=False)
ax.set_xticks([], minor=True)
ax.legend()
fig.savefig(os.path.join(f["figures"], "check_1_marginals.png"), dpi=150, bbox_inches="tight")
ratio = pd.DataFrame({f"ratio_{T}": (amax_s >= T).mean() / expected(T) for T in ck["ratio_T"]})
ratio.insert(0, "name", sites["name"])
print("\nCheck 1: simulated / expected share of years with annual max RP >= T; expected years in "
      f"{len(sim_years)}: " + ", ".join(f"T={T}: {expected(T) * len(sim_years):.0f}" for T in ck["ratio_T"]))
print(ratio.sort_values(f"ratio_{ck['ratio_T'][1]}").round(2).to_string())
for label, amax in (("simulated", amax_s), ("historical record (reference)", amax_h)):
    print(f"check 1 summary, {label}: median (range) / pooled ratio " + ", ".join(
        f"T={T} {r.median():.2f} ({r.min():.2f}-{r.max():.2f}) / {(amax >= T).values.sum() / (amax.size * expected(T)):.2f}"
        for T in ck["ratio_T"] for r in [(amax >= T).mean() / expected(T)]))

# Check 2: co-occurrence rate against distance.
lon, lat = sites["lon"].values, sites["lat"].values
dist = km(lon[:, None], lat[:, None], lon[None], lat[None])
off = ~np.eye(len(site_ids), dtype=bool)
fig, axes = plt.subplots(1, len(ck["cor_T"]), figsize=(6 * len(ck["cor_T"]), 4.5), sharey=True)
print("\nCheck 2: co-occurrence rate COR(T), simulated - historical")
for ax, T in zip(axes, ck["cor_T"]):
    cor = {}
    for label, ev, colour in (("historical", ev_h, HIST), ("simulated", ev_s, SIM)):
        e = (ev.values >= T).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            cor[label] = (e.T @ e) / e.sum(axis=0)[:, None]
        ax.scatter(dist[off], cor[label][off], s=4, alpha=0.3, color=colour, label=label)
    ok = off & np.isfinite(cor["historical"])
    diff = cor["simulated"] - cor["historical"]
    big = ok & (cor["simulated"] > ck["cor_min"])
    print(f"  T={T}: CORE {np.sqrt(np.mean(diff[ok] ** 2)):.3f} over {ok.sum()} pairs "
          f"({(off & ~ok).sum()} undefined historically), {np.sqrt(np.mean(diff[big] ** 2)):.3f} over {big.sum()} "
          f"pairs with simulated COR > {ck['cor_min']}; mean difference {diff[ok].mean():+.3f}")
    ax.set(xlabel="distance between sites (km)", title=f"COR({T})")
axes[0].set_ylabel("share of events with RP >= T at i that also have RP >= T at j")
axes[0].legend(markerscale=3)
fig.suptitle("Check 2: event co-occurrence rate")
fig.savefig(os.path.join(f["figures"], "check_2_cooccurrence.png"), dpi=150, bbox_inches="tight")

# Check 3: event extent per event and per year.
k = np.arange(len(site_ids) + 1)
fig, axes = plt.subplots(2, len(ck["extent_T"]), figsize=(4.5 * len(ck["extent_T"]), 8), sharey="row",
                         layout="constrained")
print("\nCheck 3: number of sites with RP >= T (mean / 95th percentile / max), historical vs simulated")
for j, T in enumerate(ck["extent_T"]):
    for row, (unit, h, s) in enumerate((("event", ev_h, ev_s), ("year", amax_h, amax_s))):
        n = {label: (x.values >= T).sum(axis=1) for label, x in (("historical", h), ("simulated", s))}
        for label, colour in (("historical", HIST), ("simulated", SIM)):
            axes[row, j].step(k, [(n[label] >= i).mean() for i in k], where="post", color=colour, label=label)
        axes[row, j].set(yscale="log", xlabel=f"sites with RP >= {T}", title=f"per {unit}, T = {T}")
        print(f"  per {unit}, T={T}: " + "; ".join(
            f"{label} {v.mean():.1f} / {np.percentile(v, 95):.0f} / {v.max()}" for label, v in n.items()))
axes[0, 0].set_ylabel("share of events with at least this many sites")
axes[1, 0].set_ylabel("share of years with at least this many sites")
axes[0, 0].legend()
fig.suptitle("Check 3: event extent")
fig.savefig(os.path.join(f["figures"], "check_3_extent.png"), dpi=150, bbox_inches="tight")

# Check 4: big years in the simulated annual distribution.
table = pd.DataFrame({"year": ck["big_years"], "n_events": n_hist[np.searchsorted(years, ck["big_years"])]})
fig, axes = plt.subplots(1, len(ck["extent_T"]), figsize=(4.5 * len(ck["extent_T"]), 4))
for ax, T in zip(axes, ck["extent_T"]):
    sim_n = (amax_s.values >= T).sum(axis=1)
    table[f"sites_rp_ge_{T}"] = [(amax_h.loc[y] >= T).sum() for y in ck["big_years"]]
    table[f"sim_share_ge_{T}"] = [(sim_n >= n).mean() for n in table[f"sites_rp_ge_{T}"]]
    ax.hist(sim_n, bins=np.arange(len(site_ids) + 2) - 0.5, density=True, color=SIM, alpha=0.7)
    for n, ys in table.groupby(f"sites_rp_ge_{T}")["year"]:
        ax.axvline(n, color=HIST, lw=1)
        ax.text(n, ax.get_ylim()[1], ", ".join(map(str, ys)), rotation=90, va="top", ha="right", fontsize=8)
    ax.set(xlabel=f"sites with annual max RP >= {T}", ylabel="share of simulated years", title=f"T = {T}")
fig.suptitle("Check 4: big years against the simulated annual distribution")
fig.savefig(os.path.join(f["figures"], "check_4_big_years.png"), dpi=150, bbox_inches="tight")
print("\nCheck 4: big years; sim_share_ge_T = share of simulated years with at least as many sites at RP >= T")
print(table.round(3).to_string(index=False))

# Check 5: GloFAS thresholds.
if f["rp_thresholds"] is not None:
    raise ValueError("files.rp_thresholds is set but check 5 is not written yet")
print("\nCheck 5 skipped: files.rp_thresholds is null")

# Check 6: profile of one simulated event.
e = ck["profile_event"]
info = events_sim.set_index("event_id").loc[e]
if info["source"] != "model":
    raise ValueError(f"event {e} is {info['source']}, so it has no conditioning site")
r = ev_s.loc[e]
raw = pd.read_parquet(f["sim_probs"])
cond = raw.loc[raw["event_id"] == e].set_index("site_id")["p_event"].idxmax()
dep = pd.read_parquet(f["dependence_params"])
dep = dep[dep["cond_site"] == cond].set_index("site_id").reindex(site_ids[site_ids != cond])
if dep["a"].isna().any():
    raise ValueError(f"dependence_params missing sites for conditioning site {cond}")
dep["distance_km"] = km(lon[site_ids == cond], lat[site_ids == cond], sites.loc[dep.index, "lon"], sites.loc[dep.index, "lat"])
dep["rp_years"] = r[dep.index]
dep["name"] = sites.loc[dep.index, "name"]
print(f"\nCheck 6: event {e} (simulated year {info['sim_year']}): conditioning site {cond} "
      f"{sites.loc[cond, 'name'] or '(unnamed)'} at {lon[site_ids == cond][0]}, {lat[site_ids == cond][0]}, "
      f"upstream area {sites.loc[cond, 'upstream_area_km2']:.0f} km2, rp_years {r[cond]:.0f}; "
      f"{info['n_sites_rp_ge_2']} sites with RP >= 2, {info['n_sites_rp_ge_10']} with RP >= 10")
print(f"alpha from site {cond}: median {dep['a'].median():.2f}, {(dep['a'] > 0.5).sum()} of {len(dep)} sites above 0.5; "
      f"Spearman correlation of alpha with rp_years in this event {dep['a'].corr(dep['rp_years'], method='spearman'):.2f}")
print(dep[dep["name"] != ""][["name", "distance_km", "a", "b", "rp_years"]].round(2).to_string())

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6.5))
for polygon in json.load(open(f["boundary"]))["features"][0]["geometry"]["coordinates"]:
    ax1.plot(*np.array(polygon[0]).T, color="0.5", lw=0.8)
norm = LogNorm(1, r.max())
sc = ax1.scatter(lon, lat, c=r.values, cmap="YlOrRd", norm=norm, s=40, edgecolors="0.3", linewidths=0.5)
ax1.scatter(lon[site_ids == cond], lat[site_ids == cond], marker="*", s=300, facecolors="none", edgecolors="black",
            label=f"conditioning site {cond}")
ax1.set(aspect="equal", xlabel="longitude", ylabel="latitude", title=f"rp_years in simulated event {e}")
ax1.legend(loc="lower left")
fig.colorbar(sc, ax=ax1, label="rp_years")
ax2.scatter(dep["distance_km"], dep["a"], c=dep["rp_years"], cmap="YlOrRd", norm=norm, s=40, edgecolors="0.3",
            linewidths=0.5, label="alpha (colour = rp_years in event)")
ax2.scatter(dep["distance_km"], dep["b"], marker="x", s=20, color="0.5", label="beta")
for i, row in dep[dep["name"] != ""].iterrows():
    ax2.annotate(row["name"], (row["distance_km"], row["a"]), fontsize=8, xytext=(4, 4), textcoords="offset points")
ax2.set(xlabel=f"distance from site {cond} (km)", ylabel="Heffernan-Tawn parameter",
        title=f"Dependence on site {cond}")
ax2.legend(loc="lower left")
fig.savefig(os.path.join(f["figures"], "check_6_event_profile.png"), dpi=150, bbox_inches="tight")

# Check 7: repeat flooding within a year.
T = ck["repeat_T"]
print(f"\nCheck 7: site-years with RP >= {T} in two or more events")
for label, d, index in (("historical", hist, years), ("simulated", sim, sim_years)):
    n = d[d["rp_years"] >= T].groupby(["year", "site_id"]).size().unstack().reindex(index=index, columns=site_ids).fillna(0)
    per_site = (n >= 2).sum() / (n >= 1).sum()
    print(f"  {label}: {(n >= 2).values.mean():.3f} of all site-years, {(n >= 2).values.sum() / (n >= 1).values.sum():.3f} "
          f"of site-years with at least one; per site {per_site.median():.2f} ({per_site.min():.2f}-{per_site.max():.2f})")
print(f"\nwrote check figures to {f['figures']}")
