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
import os

import cdsapi
import yaml

cfg = yaml.safe_load(open("config.yaml"))
g = cfg["glofas"]

client = cdsapi.Client()

for year in range(g["start_year"], g["end_year"] + 1):
    path = cfg["files"]["raw_discharge"].format(year=year)
    if os.path.exists(path):
        print(f"{year}: exists, skipping")
        continue
    os.makedirs(os.path.dirname(path), exist_ok=True)

    # Keys and values copied from the EWDS form's "API request" box.
    request = {
        "system_version": [g["system_version"]],
        "hydrological_model": ["lisflood"],
        "product_type": ["consolidated"],
        "timespan": ["time_mean"],
        "variable": ["average_river_discharge_in_the_last_24_hours"],
        "year": [str(year)],
        "month": [f"{m:02d}" for m in range(1, 13)],
        "day": [f"{d:02d}" for d in range(1, 32)],
        "area": g["area"],
        "data_format": "netcdf",
        "download_format": "unarchived",
    }

    # Download to a temporary name so an interrupted run never leaves a file that gets skipped.
    tmp = path + ".part"
    print(f"{year}: requesting")
    client.retrieve("cems-glofas-historical", request).download(tmp)

    # "unarchived" still returns a zip if EWDS splits the netCDF conversion into several files.
    with open(tmp, "rb") as f:
        magic = f.read(8)
    if magic != b"\x89HDF\r\n\x1a\n":
        raise ValueError(f"{year}: download is not a single netCDF-4 file (starts {magic!r}); see {tmp}")
    os.replace(tmp, path)
    print(f"{year}: wrote {path}")
