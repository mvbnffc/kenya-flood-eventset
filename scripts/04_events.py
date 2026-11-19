"""
04_events.py: identify historical flood events and each site's peak in each event.

Reads:  config.yaml (events.*, qc.max_zero_share, qc.min_mean_flow_m3s), files.discharge,
        files.site_qc, files.sites
Writes: files.events (event_id, start, end, n_days, n_sites_triggered, season),
        files.event_peaks (rows = event_id, columns = site_id, value = peak discharge)

Method:
- Use only sites with daily zero_share <= qc.max_zero_share and mean daily flow
  >= qc.min_mean_flow_m3s; print the others. Stop if a required site is dropped.
- Trigger threshold per site: the larger of (a) the daily value at which the site
  averages events.trigger_clusters_per_year exceedance clusters per year (clusters
  separated by at least events.merge_gap_days), found by stepping down through the
  site's distinct daily values from the top, the first that reaches the rate; and
  (b) the site's events.trigger_quantile_floor daily quantile. Where no value
  reaches the rate, (b) is used. Print how many sites take (b).
- A trigger day is any day when at least one site is at or above its trigger.
- Merge runs of trigger days separated by fewer than events.merge_gap_days.
  Pad each run by events.pad_days on both sides and merge runs that then overlap.
- Windows longer than events.split_above_days are split at valleys if any exist,
  otherwise kept. A valley is an interior day with fewer sites above trigger than
  some earlier day and some later day in the window. Split at the valley with the
  fewest sites above trigger, ties broken by the lowest sum over sites of
  discharge / trigger; the valley day starts the second piece. Repeat on each
  piece longer than events.split_above_days. Fail if a piece has no trigger day.
  Print every window still longer than events.split_above_days.
- Event peak at each site = maximum discharge in the window, for all sites,
  whether or not they were triggered.
- season = MAM, OND or other, from the day in the window with the most sites
  above trigger.
- Print: events per year (lambda), the window-length distribution, the split
  windows, the ten events with the most sites triggered (1997-98, 2019-20 and
  2023-24 should be among them), and the events in events.check_periods with
  their peaks at the required sites and events.check_sites.
- For each pair in events.lag_check_pairs, print the lag (days) of peak
  cross-correlation; pad_days should be at least the largest lag.
"""
import os

import numpy as np
import pandas as pd
import yaml

cfg = yaml.safe_load(open("config.yaml"))
ev, f = cfg["events"], cfg["files"]
gap, pad, split_above = ev["merge_gap_days"], ev["pad_days"], ev["split_above_days"]


def runs(above):
    """First and last positions of runs of True, merging runs separated by fewer than merge_gap_days False days."""
    idx = np.flatnonzero(above)
    breaks = np.flatnonzero(np.diff(idx) > gap)
    return idx[np.r_[0, breaks + 1]], idx[np.r_[breaks, len(idx) - 1]]


dis = pd.read_parquet(f["discharge"])
sites = pd.read_parquet(f["sites"]).merge(pd.read_parquet(f["site_qc"]), on="site_id")
cols = ["site_id", "name", "lon", "lat", "upstream_area_km2", "zero_share", "mean_flow_m3s"]
dry = sites["zero_share"] > cfg["qc"]["max_zero_share"]
low = ~dry & (sites["mean_flow_m3s"] < cfg["qc"]["min_mean_flow_m3s"])
print(f"{dry.sum()} of {len(sites)} sites excluded (zero_share > {cfg['qc']['max_zero_share']}):")
print(sites[dry][cols].round(3).to_string(index=False))
print(f"\n{low.sum()} more excluded (mean flow < {cfg['qc']['min_mean_flow_m3s']} m3/s):")
print(sites[low][cols].round(3).to_string(index=False))
if (sites["name"][dry | low] != "").any():
    raise ValueError(f"required sites excluded: {list(sites['name'][(dry | low) & (sites['name'] != '')])}")
sites = sites[~(dry | low)]
dis = dis[sites["site_id"].astype(str)]
q = dis.values
n_years = dis.index.year.nunique()

# Trigger threshold per site.
target = ev["trigger_clusters_per_year"] * n_years
floor = np.quantile(q, ev["trigger_quantile_floor"], axis=0)
trigger = floor.copy()
for j in range(q.shape[1]):
    for u in np.unique(q[:, j])[::-1]:
        if u <= floor[j]:
            break
        if len(runs(q[:, j] >= u)[0]) >= target:
            trigger[j] = u
            break
print(f"\n{(trigger == floor).sum()} of {len(trigger)} sites use the trigger_quantile_floor "
      f"({ev['trigger_quantile_floor']}) quantile")
if (trigger <= 0).any():
    raise ValueError(f"zero trigger at sites {list(dis.columns[trigger <= 0])}")
quantile = (q < trigger).mean(axis=0)
print(f"\n{len(trigger)} sites used; trigger quantile min {quantile.min():.4f}, "
      f"median {np.median(quantile):.4f}, max {quantile.max():.4f}")

exceed = q >= trigger
n_above = exceed.sum(axis=1)
trig_day = n_above > 0
relative = (q / trigger).sum(axis=1)

# Windows: merged runs of trigger days, padded, overlapping windows merged.
windows = []
for a, b in zip(*runs(trig_day)):
    a, b = max(a - pad, 0), min(b + pad, len(q) - 1)
    if windows and a <= windows[-1][1]:
        windows[-1][1] = b
    else:
        windows.append([a, b])

# Split long windows at valleys, repeatedly; a long piece with no valley is kept.
pieces, split = [], []
for a, b in windows:
    todo, done = [(a, b)], []
    while todo:
        a2, b2 = todo.pop()
        n = n_above[a2:b2 + 1]
        interior = np.arange(a2 + 1, b2)
        before = np.maximum.accumulate(n)[interior - a2 - 1]
        after = np.maximum.accumulate(n[::-1])[::-1][interior - a2 + 1]
        valleys = interior[(n_above[interior] < before) & (n_above[interior] < after)]
        if b2 - a2 + 1 <= split_above or len(valleys) == 0:
            done.append((a2, b2))
            continue
        t = valleys[np.lexsort((relative[valleys], n_above[valleys]))[0]]
        todo += [(a2, t - 1), (t, b2)]
    done.sort()
    if len(done) > 1:
        split.append((b - a + 1, [int(d2 - d1 + 1) for d1, d2 in done]))
    pieces += done
for a, b in pieces:
    if not trig_day[a:b + 1].any():
        raise ValueError(f"window {dis.index[a].date()} to {dis.index[b].date()} has no trigger day")

dates = dis.index
busiest = [a + n_above[a:b + 1].argmax() for a, b in pieces]
months = dates[busiest].month
events = pd.DataFrame({
    "event_id": np.arange(1, len(pieces) + 1),
    "start": dates[[a for a, _ in pieces]],
    "end": dates[[b for _, b in pieces]],
    "n_days": [b - a + 1 for a, b in pieces],
    "n_sites_triggered": [exceed[a:b + 1].any(axis=0).sum() for a, b in pieces],
    "season": np.where(months.isin([3, 4, 5]), "MAM", np.where(months.isin([10, 11, 12]), "OND", "other")),
})
peaks = pd.DataFrame([q[a:b + 1].max(axis=0) for a, b in pieces], columns=dis.columns,
                     index=pd.Index(events["event_id"], name="event_id"))

os.makedirs(os.path.dirname(f["events"]), exist_ok=True)
events.to_parquet(f["events"], index=False)
peaks.to_parquet(f["event_peaks"])
print(f"wrote {f['events']} and {f['event_peaks']}: {len(events)} events x {peaks.shape[1]} sites")

pd.set_option("display.width", 200)
print(f"\nlambda = {len(events)} events / {n_years} years = {len(events) / n_years:.2f} per year")
print("season counts:", events["season"].value_counts().to_dict())
print("window length (days):", events["n_days"].describe().round(1).to_dict())
print(f"\n{len(split)} of {len(windows)} windows split; the five longest (before -> after):")
for before, after in sorted(split, reverse=True)[:5]:
    print(f"  {before} -> {after}")
long = events[events["n_days"] > split_above]
print(f"\n{len(long)} windows still longer than {split_above} days (no valley):")
print(long[["event_id", "start", "end", "n_days", "n_sites_triggered"]].to_string(index=False))

print("\nTen events with the most sites triggered:")
print(events.sort_values(["n_sites_triggered", "start"], ascending=[False, True]).head(10).to_string(index=False))

shown = sites[(sites["name"] != "") | sites["site_id"].isin(ev["check_sites"])]
labels = {str(i): n or str(i) for i, n in zip(shown["site_id"], shown["name"])}
for label, first, last in ev["check_periods"]:
    sel = (events["start"] <= pd.Timestamp(last)) & (events["end"] >= pd.Timestamp(first))
    table = events[sel][["event_id", "start", "end", "n_sites_triggered"]].merge(
        peaks[list(labels)].rename(columns=labels).round(0), left_on="event_id", right_index=True)
    print(f"\n{label} ({first} to {last}): peak discharge (m3/s)")
    print(table.to_string(index=False) if len(table) else "  no events")

print(f"\nLag of peak cross-correlation (pad_days = {pad}):")
for up, down in ev["lag_check_pairs"]:
    x, y = dis[str(up)], dis[str(down)]
    corr = [y.corr(x.shift(lag)) for lag in range(split_above + 1)]
    print(f"  {up} -> {down}: {int(np.argmax(corr))} days (r = {max(corr):.3f})")
