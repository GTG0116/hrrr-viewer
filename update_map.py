import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from herbie import Herbie
from datetime import datetime, timedelta
import pytz
import os
import numpy as np
import xarray as xr
from matplotlib.colors import ListedColormap

# --- CONFIGURATION ---
ZOOM_BOUNDS = [-82, -67, 37, 48] # [West, East, South, North]

# 1. Find the most recent available HRRR run
H_init = None
for offset in range(24):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        date_str = search_time.strftime('%Y-%m-%d %H:00')
        print(f"Checking for data at: {date_str}...")
        
        test_H = Herbie(date_str, model='hrrr', product='sfc', fxx=0, priority=['aws'])
        if test_H.inventory() is not None:
            H_init = test_H
            print(f"Success! Using run: {date_str}")
            break
    except:
        continue

if not H_init:
    raise Exception("Could not find any valid HRRR data on AWS.")

os.makedirs("frames_temp", exist_ok=True)
os.makedirs("frames_precip", exist_ok=True)

init_hour = int(H_init.date.strftime('%H'))
max_fxx = 48 if init_hour in [0, 6, 12, 18] else 18

# 2. Loop through frames
for fxx in range(max_fxx + 1):
    try:
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), 
                   model='hrrr', product='sfc', fxx=fxx, priority=['aws'])
        
        # Pull data and force it into a single Dataset
        ds = H.xarray("(:TMP:2 m|:CSNOW:|:CICEP:|:CFRZR:|:CRAIN:)")
        if isinstance(ds, list):
            ds = xr.merge(ds)

        # Get the projection from Herbie directly, which is more reliable
        map_projection = H.crs 

        utc_v = H.valid_date.replace(tzinfo=pytz.UTC)
        et_v = utc_v.astimezone(pytz.timezone('US/Eastern'))
        timestamp = f"Valid: {utc_v.strftime('%m/%d %H:%M')}Z | {et_v.strftime('%I:%M %p ET')}"

        # --- MAP 1: TEMPERATURE ---
        temp_f = (ds.t2m - 273.15) * 9/5 + 32
        fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection': map_projection})
        ax.set_extent(ZOOM_BOUNDS, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        
        im = ax.pcolormesh(ds.longitude, ds.latitude, temp_f, transform=ccrs.PlateCarree(), 
                          cmap='jet', vmin=0, vmax=100)
        plt.colorbar(im, label="Temperature (°F)", orientation='horizontal', pad=0.05)
        plt.title(f"HRRR Temperature f{fxx:02d}\n{timestamp}", loc='left', fontweight='bold')
        plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()

        # --- MAP 2: PRECIP TYPE ---
        precip_mask = np.zeros_like(temp_f)
        if 'crain' in ds: precip_mask = np.where(ds.crain == 1, 1, precip_mask)
        if 'csnow' in ds: precip_mask = np.where(ds.csnow == 1, 2, precip_mask)
        if 'cfrzr' in ds: precip_mask = np.where(ds.cfrzr == 1, 3, precip_mask)
        if 'cicep' in ds: precip_mask = np.where(ds.cicep == 1, 4, precip_mask)

        fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection': map_projection})
        ax.set_extent(ZOOM_BOUNDS, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        
        p_cmap = ListedColormap(['none', '#2ecc71', '#3498db', '#e74c3c', '#e67e22'])
        ax.pcolormesh(ds.longitude, ds.latitude, np.ma.masked_where(precip_mask == 0, precip_mask),
                     transform=ccrs.PlateCarree(), cmap=p_cmap, vmin=0, vmax=4)
        
        plt.title(f"HRRR Precip Type f{fxx:02d}\n{timestamp}", loc='left', fontweight='bold')
        plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()
        
        print(f"Success: f{fxx:02d}")

    except Exception as e:
        print(f"Skipping f{fxx:02d}: {e}")
