"""
04_events.py: identify historical flood events and each site's peak in each event.

Reads:  config.yaml (events.*), files.discharge
Writes: files.events (event_id, start, end, n_days, n_sites_triggered, season),
        files.event_peaks (rows = event_id, columns = site_id, value = peak discharge)

Method:
- Trigger threshold per site: the daily quantile at which the site averages
  events.trigger_clusters_per_year exceedance clusters per year (clusters separated
  by at least events.merge_gap_days). Find it by a simple search over quantiles.
- A trigger day is any day when at least one site is above its trigger.
- Merge runs of trigger days separated by fewer than events.merge_gap_days.
  Pad each run by events.pad_days on both sides and merge runs that then overlap.
  Split any window longer than events.max_length_days into consecutive windows.
- Event peak at each site = maximum discharge in the window, for all sites,
  whether or not they were triggered.
- season = MAM, OND or other, from the day in the window with the most sites
  above trigger.
- Print: events per year (lambda), the window-length distribution, and the ten
  events with the most sites triggered (1997-98, 2019-20 and 2023-24 should be
  among them).
- For each pair in events.lag_check_pairs, print the lag (days) of peak
  cross-correlation; pad_days should be at least the largest lag.
"""
import yaml

cfg = yaml.safe_load(open("config.yaml"))

# To be written in an agent session.
