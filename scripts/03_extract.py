"""
03_extract.py: extract daily discharge at each site, and basic QC.

Reads:  config.yaml, files.raw_discharge (all years), files.sites
Writes: files.discharge (rows = dates, columns = site_id),
        files.site_qc (one row per site)

Method:
- Open one year at a time, select the site cells by row and col, and concatenate.
  Fail if any year is missing or the grid differs between years.
- Site QC columns: site_id, n_missing, n_negative, zero_share (share of days with
  discharge < 0.01 m3/s), annual_max_trend (Spearman correlation of annual maximum
  with year) and annual_max_trend_p.
- Print the sites with missing or negative values, or a trend with p < 0.05.
"""
import yaml

cfg = yaml.safe_load(open("config.yaml"))

# To be written in an agent session.
