"""
01_download.py: download GloFAS v5 historical daily river discharge for Kenya.

Reads:  config.yaml (glofas.*, files.raw_discharge)
Writes: files.raw_discharge, one netCDF file per year

Method:
- One cdsapi request per year from glofas.start_year to glofas.end_year,
  area = glofas.area, system version = glofas.system_version, consolidated data only.
- Copy the request keys from the "API request" box on the EWDS download form
  rather than guessing them.
- Request netCDF. If only GRIB is offered, stop and ask before adding cfgrib.
- Skip years whose file already exists, so the script can be re-run after a dropout.
"""
import yaml

cfg = yaml.safe_load(open("config.yaml"))

# To be written in an agent session.
