"""
02_sites.py: choose the river sites (GloFAS cells) used as virtual gauges.

Reads:  config.yaml (sites.*), files.upstream_area, files.ldd, files.lake_mask,
        files.boundary, files.regulated_sites,
        files.raw_discharge (first year, grid check only)
Writes: files.sites, figures/sites.png

Method:
- Candidate cells: upstream area >= sites.min_upstream_area_km2, inside Kenya,
  not a lake cell (lake mask == 1).
- Use the flow-direction (LDD) map to follow each river downstream. Place a site
  at each confluence of two candidate rivers, wherever upstream area has grown by
  sites.area_growth_for_new_site since the last site, and at each river's last
  candidate cell: no candidate cell anywhere further downstream. A river that
  crosses the border several times gets one exit site, where it leaves Kenya for
  the last time; a river entering a lake ends at the shore.
- Drop every site within sites.excluded_radius_cells of a point in
  sites.excluded_sites (e.g. the Omo delta, whose flow is from Ethiopia).
- Add sites.required_sites, snapped to the nearest candidate cell.
- If sites.max_sites is set, keep all required sites, then fill up to max_sites
  with the largest upstream area (for the small end-to-end run).
- Columns: site_id, name (blank unless required), lon, lat, row, col,
  upstream_area_km2, outlet_id (basin label), regulated (downstream of a
  reservoir outlet in files.regulated_sites).
- Print upstream area at the required sites so they can be checked by hand
  against published WRA catchment areas.
- Plot sites on the river network, sized by upstream area.
"""
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from matplotlib.path import Path

cfg = yaml.safe_load(open("config.yaml"))
s, f = cfg["sites"], cfg["files"]
N, W, S, E = cfg["glofas"]["area"]

dis = xr.open_dataset(f["raw_discharge"].format(year=cfg["glofas"]["start_year"]))
lat, lon = dis.latitude.values, dis.longitude.values
nr, nc = len(lat), len(lon)


def load_grid(path):
    """Cut a global GloFAS auxiliary map to glofas.area; fail unless it matches the discharge grid."""
    band = xr.open_dataset(path).Band1.sel(lat=slice(N, S), lon=slice(W, E))
    if band.shape != (nr, nc) or not (np.allclose(band.lat, lat, atol=1e-6) and np.allclose(band.lon, lon, atol=1e-6)):
        raise ValueError(f"{path}: grid does not match the discharge grid over glofas.area")
    return band.values


def downstream(i):
    """Flat indices of cell i and every cell below it, to a pit or the edge of glofas.area."""
    path = [i]
    while down[path[-1]] >= 0:
        path.append(down[path[-1]])
    return path


def cells_from(x, y):
    """Distance in grid cells from the point (lon x, lat y) to every cell centre."""
    return np.hypot(centres[:, 0] - x, centres[:, 1] - y) / abs(lon[1] - lon[0])


up =load_grid(f["upstream_area"]).ravel() / 1e6  # m2 -> km2
ldd = load_grid(f["ldd"])
lake_values = load_grid(f["lake_mask"]).ravel()
if not set(np.unique(lake_values[~np.isnan(lake_values)])) <= {0, 1}:
    raise ValueError(f"{f['lake_mask']}: expected values 0 and 1")
lake = lake_values == 1
reservoirs = load_grid(f["regulated_sites"]).ravel()

# Downstream neighbour of each cell (flat index, -1 at pits and where flow leaves the area).
# PCRaster LDD codes are numeric-keypad directions; 5 is a pit. Rows run north to south.
codes = set(np.unique(ldd[~np.isnan(ldd)]).astype(int))
if not codes <= set(range(1, 10)):
    raise ValueError(f"unexpected LDD codes: {sorted(codes)}")
step = {1: (1, -1), 2: (1, 0), 3: (1, 1), 4: (0, -1), 6: (0, 1), 7: (-1, -1), 8: (-1, 0), 9: (-1, 1)}
rows, cols = np.indices((nr, nc))
dr, dc = np.zeros((nr, nc), int), np.zeros((nr, nc), int)
for code, (a, b) in step.items():
    dr[ldd == code], dc[ldd == code] = a, b
rr, cc = rows + dr, cols + dc
flows = ((dr != 0) | (dc != 0)) & (rr >= 0) & (rr < nr) & (cc >= 0) & (cc < nc)
down = np.where(flows, rr * nc + cc, -1).ravel()

# Kenya mask from cell centres.
features = json.load(open(f["boundary"]))["features"]
if len(features) != 1 or features[0]["geometry"]["type"] != "MultiPolygon":
    raise ValueError(f"{f['boundary']}: expected one MultiPolygon feature")
lon2, lat2 = np.meshgrid(lon, lat)
centres = np.column_stack([lon2.ravel(), lat2.ravel()])
kenya = np.zeros(nr * nc, bool)
for polygon in features[0]["geometry"]["coordinates"]:
    if len(polygon) != 1:
        raise ValueError(f"{f['boundary']}: polygon with holes")
    kenya |= Path(polygon[0]).contains_points(centres)

cand = kenya & ~lake & (up >= s["min_upstream_area_km2"])
has_down = down >= 0
n_in = np.bincount(down[cand & has_down], minlength=nr * nc)  # candidate cells flowing in
down_cand = np.zeros(nr * nc, bool)
down_cand[has_down] = cand[down[has_down]]
confluence = cand & (n_in >= 2)
last = cand & ~down_cand
crossings = last.sum()
for i in np.flatnonzero(last):
    last[i] = not cand[downstream(i)[1:]].any()
print(f"{cand.sum()} candidate cells, {confluence.sum()} confluences, {last.sum()} last cells "
      f"({crossings - last.sum()} earlier border crossings dropped)")

# Walk downstream. ref = upstream area at the last site above the cell; a river's
# first candidate cell is its starting reference. Sorting by upstream area visits
# every cell after all the cells that drain into it.
site = confluence | last
ref = up.copy()
growth = 1 + s["area_growth_for_new_site"]
for i in np.flatnonzero(cand)[np.argsort(up[cand])]:
    if up[i] >= growth * ref[i]:
        site[i] = True
    j = down[i]
    if j >= 0 and cand[j] and n_in[j] == 1:
        ref[j] = up[i] if site[i] else ref[i]
print(f"{site.sum()} sites from the river network")

# Excluded points: drop every site within sites.excluded_radius_cells of each.
excluded = np.zeros(nr * nc, bool)
for name, x, y in s["excluded_sites"]:
    near = site & (cells_from(x, y) <= s["excluded_radius_cells"])
    if not near.any():
        raise ValueError(f"excluded site {name}: no site within {s['excluded_radius_cells']} cells")
    print(f"excluded {name}: " + ", ".join(f"{centres[i, 1]:.3f} {centres[i, 0]:.3f} ({up[i]:.0f} km2)"
                                           for i in np.flatnonzero(near)))
    excluded |= near
site &= ~excluded

# Required sites, snapped to the nearest candidate cell.
names = {}
print("\nRequired sites (check upstream area against WRA catchment areas):")
print(f"  {'name':16s} {'lat':>7s} {'lon':>7s}  {'upstream':>13s}  lake  snapped")
for name, x, y in s["required_sites"]:
    dist = np.where(cand, cells_from(x, y), np.inf)
    i = int(dist.argmin())
    if i in names:
        raise ValueError(f"{name} and {names[i]} snap to the same cell")
    if excluded[i]:
        raise ValueError(f"{name} snaps to an excluded site")
    names[i] = name
    site[i] = True
    print(f"  {name:16s} {centres[i, 1]:7.3f} {centres[i, 0]:7.3f}  {up[i]:9.0f} km2  {str(lake[i]):5s} "
          f"{dist[i] * abs(lon[1] - lon[0]) * 111:.1f} km")

# site_id ranks all sites by upstream area, so ids do not change with max_sites.
idx = np.flatnonzero(site)
idx = idx[np.argsort(-up[idx], kind="stable")]
sites = pd.DataFrame({
    "site_id": np.arange(1, len(idx) + 1),
    "name": [names.get(i, "") for i in idx],
    "lon": np.round(centres[idx, 0], 4),
    "lat": np.round(centres[idx, 1], 4),
    "row": idx // nc,
    "col": idx % nc,
    "upstream_area_km2": up[idx],
    "outlet_id": [downstream(i)[-1] for i in idx],
})

regulated = set()
print("\nReservoir outlets (files.regulated_sites):")
for i in np.flatnonzero(~np.isnan(reservoirs)):
    print(f"  id {reservoirs[i]:.0f} at {centres[i, 1]:.3f} {centres[i, 0]:.3f}, {up[i]:.0f} km2")
    regulated.update(downstream(i))
sites["regulated"] = [i in regulated for i in idx]

if s["max_sites"] is not None:
    required = sites[sites["name"] != ""]
    if len(required) > s["max_sites"]:
        raise ValueError(f"{len(required)} required sites exceed max_sites = {s['max_sites']}")
    others = sites[sites["name"] == ""].head(s["max_sites"] - len(required))
    sites = pd.concat([required, others]).sort_values("site_id").reset_index(drop=True)

print(f"\n{len(sites)} sites kept, {sites['regulated'].sum()} regulated, "
      f"{sites['outlet_id'].nunique()} outlets")
print(sites[sites["regulated"]][["site_id", "name", "lon", "lat", "upstream_area_km2"]].to_string(index=False))

os.makedirs(os.path.dirname(f["sites"]), exist_ok=True)
sites.to_parquet(f["sites"], index=False)
print(f"wrote {f['sites']}")

# Figure: network shaded by upstream area, sites sized by upstream area.
fig, ax = plt.subplots(figsize=(8, 9))
network = np.where(up >= s["min_upstream_area_km2"], np.log10(up), np.nan).reshape(nr, nc)
extent = [lon[0] - 0.025, lon[-1] + 0.025, lat[-1] - 0.025, lat[0] + 0.025]
ax.imshow(np.where(lake, 1, np.nan).reshape(nr, nc), extent=extent, cmap="Greys", vmin=0, vmax=4,
          interpolation="nearest")
im = ax.imshow(np.where(lake.reshape(nr, nc), np.nan, network), extent=extent, cmap="Blues",
               interpolation="nearest")
for polygon in features[0]["geometry"]["coordinates"]:
    ring = np.array(polygon[0])
    ax.plot(ring[:, 0], ring[:, 1], color="0.4", lw=0.8)
size = 4 + 60 * (sites["upstream_area_km2"] / sites["upstream_area_km2"].max()) ** 0.5
for reg, colour, label in [(False, "#1b1b1b", "site"), (True, "#d95f02", "site, regulated")]:
    sel = sites["regulated"] == reg
    ax.scatter(sites["lon"][sel], sites["lat"][sel], s=size[sel], c=colour, edgecolors="white",
               linewidths=0.5, label=label, zorder=3)
res = np.flatnonzero(~np.isnan(reservoirs))
ax.scatter(centres[res, 0], centres[res, 1], marker="^", s=40, facecolors="none", edgecolors="#d95f02",
           label="reservoir outlet", zorder=4)
ax.scatter([], [], marker="s", s=40, color=plt.cm.Greys(0.25), label="lake")
for _, r in sites[sites["name"] != ""].iterrows():
    ax.annotate(r["name"], (r["lon"], r["lat"]), xytext=(5, 3), textcoords="offset points", fontsize=7)
fig.colorbar(im, ax=ax, shrink=0.6, label="log10 upstream area (km2)")
ax.legend(loc="lower left", fontsize=8)
ax.set_aspect("equal")
ax.set_title(f"{len(sites)} sites, upstream area >= {s['min_upstream_area_km2']} km2")
os.makedirs(f["figures"], exist_ok=True)
fig.savefig(os.path.join(f["figures"], "sites.png"), dpi=150, bbox_inches="tight")
print(f"wrote {os.path.join(f['figures'], 'sites.png')}")
