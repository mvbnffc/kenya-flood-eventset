"""
03_extract.py: extract daily discharge at each site, and basic QC.

Reads:  config.yaml, files.raw_discharge (all years), files.sites
Writes: files.discharge (rows = dates, columns = site_id),
        files.site_qc (one row per site)

Method:
- Open one year at a time, select the site cells by row and col, and concatenate.
  Fail if any year is missing or the grid differs between years.
- Each value is the 24-hour mean ending at valid_time, so it is labelled with the
  previous day's date (1980-01-02 00:00 -> 1980-01-01). Fail unless the dates run
  daily without gaps from glofas.start_year-01-01 to glofas.end_year-12-31.
- Site QC columns: site_id, n_missing, n_negative, zero_share (share of days with
  discharge < qc.zero_discharge_m3s), annual_max_trend (Spearman correlation of
  calendar-year maximum with year) and annual_max_trend_p.
- Print the sites with missing or negative values, or a trend with p < qc.trend_p.
"""
import os

import numpy as np
import pandas as pd
import xarray as xr
import yaml
from scipy.stats import spearmanr

cfg = yaml.safe_load(open("config.yaml"))
g, f, qc = cfg["glofas"], cfg["files"], cfg["qc"]

sites = pd.read_parquet(f["sites"])
rows = xr.DataArray(sites["row"].values, dims="site")
cols = xr.DataArray(sites["col"].values, dims="site")

parts, grid = [], None
for year in range(g["start_year"], g["end_year"] + 1):
    path = f["raw_discharge"].format(year=year)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    ds = xr.open_dataset(path)
    if grid is None:
        grid = (ds.latitude.values, ds.longitude.values)
        # The raw files are already cut to glofas.area: the grid 02 indexed by row and col.
        if not (np.allclose(grid[0][sites["row"]], sites["lat"], atol=1e-4)
                and np.allclose(grid[1][sites["col"]], sites["lon"], atol=1e-4)):
            raise ValueError(f"{path}: site row/col do not match site lon/lat")
    elif not (np.array_equal(ds.latitude.values, grid[0]) and np.array_equal(ds.longitude.values, grid[1])):
        raise ValueError(f"{path}: grid differs from {g['start_year']}")
    q = ds["avg_dis"].isel(latitude=rows, longitude=cols).transpose("valid_time", "site").values
    dates = pd.DatetimeIndex(ds["valid_time"].values) - pd.Timedelta(days=1)
    parts.append(pd.DataFrame(q, index=dates, columns=sites["site_id"].astype(str)))
    ds.close()

dis = pd.concat(parts)
dis.index.name = "date"
expected = pd.date_range(f"{g['start_year']}-01-01", f"{g['end_year']}-12-31", freq="D")
if not dis.index.equals(expected):
    raise ValueError(f"dates are not a daily run from {expected[0].date()} to {expected[-1].date()}")

os.makedirs(os.path.dirname(f["discharge"]), exist_ok=True)
dis.to_parquet(f["discharge"])
print(f"wrote {f['discharge']}: {dis.shape[0]} rows (dates) x {dis.shape[1]} columns (sites), "
      f"{dis.index[0].date()} to {dis.index[-1].date()}")

annual_max = dis.groupby(dis.index.year).max()
trend = [spearmanr(annual_max.index, annual_max[c]) for c in dis.columns]
site_qc = pd.DataFrame({
    "site_id": sites["site_id"].values,
    "n_missing": dis.isna().sum().values,
    "n_negative": (dis < 0).sum().values,
    "zero_share": (dis < qc["zero_discharge_m3s"]).mean().values,
    "annual_max_trend": [t.statistic for t in trend],
    "annual_max_trend_p": [t.pvalue for t in trend],
})
site_qc.to_parquet(f["site_qc"], index=False)
print(f"wrote {f['site_qc']}")

report = sites[["site_id", "name", "upstream_area_km2", "regulated"]].merge(site_qc, on="site_id")
flags = {
    "missing values": report["n_missing"] > 0,
    "negative values": report["n_negative"] > 0,
    f"annual-maximum trend with p < {qc['trend_p']}": report["annual_max_trend_p"] < qc["trend_p"],
}
pd.set_option("display.width", 200)
for label, flagged in flags.items():
    print(f"\n{flagged.sum()} of {len(report)} sites with {label}")
    if flagged.any():
        print(report[flagged].round(3).to_string(index=False))

print(f"\nzero_share (days < {qc['zero_discharge_m3s']} m3/s) at the required sites:")
print(report[report["name"] != ""][["site_id", "name", "upstream_area_km2", "zero_share"]].to_string(index=False))
