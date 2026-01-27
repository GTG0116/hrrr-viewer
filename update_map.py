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
import warnings

# Suppress warnings
warnings.filterwarnings("ignore")

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
            print(f"Success! Found run: {date_str}")
            break
    except:
        continue

if not H_init:
    raise Exception("Could not find any valid HRRR data on AWS.")

# Create output folders
os.makedirs("frames_temp", exist_ok=True)
os.makedirs("frames_precip", exist_ok=True)

# Set run length
init_hour = int(H_init.date.strftime('%H'))
max_fxx = 48 if init_hour in [0, 6, 12, 18] else 18

print(f"Starting generation for {max_fxx} frames...")

# 2. Loop through frames
for fxx in range(max_fxx + 1):
    try:
        # Re-initialize Herbie
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), 
                   model='hrrr', product='sfc', fxx=fxx, priority=['aws'])

        # Time strings
        utc_v = H.valid_date.replace(tzinfo=pytz.UTC)
        et_v = utc_v.astimezone(pytz.timezone('US/Eastern'))
        timestamp = f"Valid: {utc_v.strftime('%m/%d %H:%M')}Z | {et_v.strftime('%I:%M %p ET')}"

        # --- PART A: TEMPERATURE MAP ---
        try:
            # 1. Get Temperature Data
            ds_t = H.xarray("TMP:2 m")
            if isinstance(ds_t, list): ds_t = ds_t[0]

            # 2. Extract Projection (Safely)
            try:
                map_crs = ds_t.herbie.crs
            except:
                # Fallback: Manual HRRR Projection
                map_crs = ccrs.LambertConformal(central_longitude=-97.5, 
                                                central_latitude=38.5, 
                                                standard_parallels=(38.5, 38.5))

            temp_f = (ds_t.t2m - 273.15) * 9/5 + 32

            # Plot Temp
            fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection': map_crs})
            ax.set_extent(ZOOM_BOUNDS, crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
            
            im = ax.pcolormesh(ds_t.longitude, ds_t.latitude, temp_f, transform=ccrs.PlateCarree(), 
                               cmap='jet', vmin=0, vmax=100)
            plt.colorbar(im, label="Temperature (°F)", orientation='horizontal', pad=0.05)
            plt.title(f"HRRR Temp f{fxx:02d}\n{timestamp}", loc='left', fontweight='bold')
            plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()

        except Exception as e_temp:
            print(f"Failed to create Temp map for f{fxx}: {e_temp}")
            continue # If temp fails, we probably can't do precip either, so skip to next hour

        # --- PART B: PRECIP TYPE MAP ---
        try:
            # Fetch ONLY precip types
            ds_p = H.xarray(":(CRAIN|CSNOW|CFRZR|CICEP):")
            
            precip_mask = np.zeros_like(temp_f)

            def process_precip_ds(d):
                # Ensure we are looking at the right variable keys
                mask = np.zeros_like(temp_f)
                if 'crain' in d: mask = np.where(d.crain == 1, 1, mask)
                if 'csnow' in d: mask = np.where(d.csnow == 1, 2, mask)
                if 'cfrzr' in d: mask = np.where(d.cfrzr == 1, 3, mask)
                if 'cicep' in d: mask = np.where(d.cicep == 1, 4, mask)
                return mask

            if isinstance(ds_p, list):
                for sub_ds in ds_p:
                    precip_mask = np.maximum(precip_mask, process_precip_ds(sub_ds))
            else:
                precip_mask = process_precip_ds(ds_p)

            # Plot Precip
            fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection': map_crs})
            ax.set_extent(ZOOM_BOUNDS, crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
            
            p_cmap = ListedColormap(['none', '#2ecc71', '#3498db', '#e74c3c', '#e67e22'])
            ax.pcolormesh(ds_t.longitude, ds_t.latitude, np.ma.masked_where(precip_mask == 0, precip_mask),
                          transform=ccrs.PlateCarree(), cmap=p_cmap, vmin=0, vmax=4)
            
            plt.title(f"HRRR Precip Type f{fxx:02d}\n{timestamp}", loc='left', fontweight='bold')
            plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
            
        except Exception as e_precip:
             print(f"Failed to create Precip map for f{fxx}: {e_precip}")

        print(f"Finished f{fxx:02d}")

    except Exception as e:
        print(f"Critical error on f{fxx}: {e}")
