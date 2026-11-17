"""
02_sites.py: choose the river sites (GloFAS cells) used as virtual gauges.

Reads:  config.yaml (sites.*), files.upstream_area, files.ldd
Writes: files.sites, figures/sites.png

Method:
- Candidate cells: upstream area >= sites.min_upstream_area_km2, inside Kenya.
- Use the flow-direction (LDD) map to follow each river downstream. Place a site
  at each confluence of two candidate rivers, wherever upstream area has grown by
  sites.area_growth_for_new_site since the last site, and at each river's last
  cell in Kenya.
- Add sites.required_sites, snapped to the nearest candidate cell.
- If sites.max_sites is set, keep that many sites with the largest upstream area
  (for the small end-to-end run).
- Columns: site_id, name (blank unless required), lon, lat, row, col,
  upstream_area_km2, outlet_id (basin label), regulated (from sites.regulated).
- Print upstream area at the required sites so they can be checked by hand
  against published WRA catchment areas.
- Plot sites on the river network, sized by upstream area.
"""
import yaml

cfg = yaml.safe_load(open("config.yaml"))

# To be written in an agent session.
