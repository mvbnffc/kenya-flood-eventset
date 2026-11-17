"""
06_event_rp.py: build the event-return-period lookup table.

Reads:  config.yaml (glofas.*, output.rp_floor), files.sim_probs, files.events
Writes: files.events_sim (event_id, sim_year, source, n_sites_rp_ge_2, n_sites_rp_ge_10),
        files.event_rp (event_id, site_id, p_event, aep, rp_years), only rows
        with rp_years >= output.rp_floor

Method:
- lambda = number of historical events / number of record years
  (glofas.end_year - glofas.start_year + 1).
- aep = 1 - exp(-lambda * (1 - p_event)); rp_years = 1 / aep.
- Print the number of simulated years, events and rows in event_rp.
"""
import yaml

cfg = yaml.safe_load(open("config.yaml"))

# To be written in an agent session.
