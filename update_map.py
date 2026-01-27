import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from herbie import Herbie
from datetime import datetime, timedelta
import pytz
import os
import numpy as np
from matplotlib.colors import ListedColormap

# --- CONFIGURATION ---
ZOOM_BOUNDS = [-82, -67, 37, 48] 

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
        
        # Pull Projection and Time metadata directly from Herbie object
        map_projection = H.crs 
        utc_v = H.valid_date.replace(tzinfo=pytz.UTC)
        et_v = utc_v.astimezone(pytz.timezone('US/Eastern'))
        timestamp = f"Valid: {utc_v.strftime('%m/%d %H:%M')}Z | {et_v.strftime('%I:%M %p ET')}"

        # --- STEP A: TEMPERATURE ---
        ds_t = H.xarray("TMP:2 m")
        if isinstance(ds_t, list): ds_t = ds_t[0]
        temp_f = (ds_t.t2m - 273.15) * 9/5 + 32

        fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection': map_projection})
        ax.set_extent(ZOOM_BOUNDS, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        im = ax.pcolormesh(ds_t.longitude, ds_t.latitude, temp_f, transform=ccrs.PlateCarree(), 
                          cmap='jet', vmin=0, vmax=100)
        plt.colorbar(im, label="Temperature (°F)", orientation='horizontal', pad=0.05)
        plt.title(f"HRRR Temp f{fxx:02d}\n{timestamp}", loc='left', fontweight='bold')
        plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()

        # --- STEP B: PRECIP TYPE ---
        # Fetching separately avoids the 'Multiple Hypercube' list error
        precip_mask = np.zeros_like(temp_f)
        try:
            ds_p = H.xarray(":(CRAIN|CSNOW|CFRZR|CICEP):")
            if isinstance(ds_p, list):
                # Manual merge loop for precip types
                for sub in ds_p:
                    if 'crain' in sub: precip_mask = np.where(sub.crain == 1, 1, precip_mask)
                    if 'csnow' in sub: precip_mask = np.where(sub.csnow == 1, 2, precip_mask)
                    if 'cfrzr' in sub: precip_mask = np.where(sub.cfrzr == 1, 3, precip_mask)
                    if 'cicep' in sub: precip_mask = np.where(sub.cicep == 1, 4, precip_mask)
            else:
                if 'crain' in ds_p: precip_mask = np.where(ds_p.crain == 1, 1, precip_mask)
                if 'csnow' in ds_p: precip_mask = np.where(ds_p.csnow == 1, 2, precip_mask)
                if 'cfrzr' in ds_p: precip_mask = np.where(ds_p.cfrzr == 1, 3, precip_mask)
                if 'cicep' in ds_p: precip_mask = np.where(ds_p.cicep == 1, 4, precip_mask)
        except:
            print(f"No precip data for f{fxx}")

        fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection': map_projection})
        ax.set_extent(ZOOM_BOUNDS, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        p_cmap = ListedColormap(['none', '#2ecc71', '#3498db', '#e74c3c', '#e67e22'])
        ax.pcolormesh(ds_t.longitude, ds_t.latitude, np.ma.masked_where(precip_mask == 0, precip_mask),
                     transform=ccrs.PlateCarree(), cmap=p_cmap, vmin=0, vmax=4)
        plt.title(f"HRRR Precip f{fxx:02d}\n{timestamp}", loc='left', fontweight='bold')
        plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()
        
        print(f"Success: f{fxx:02d}")

    except Exception as e:
        print(f"Skipping f{fxx:02d}: {e}")
