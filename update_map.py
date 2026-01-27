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

# Suppress warnings to keep logs clean
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
        
        # We check fxx=0 to see if the run exists
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

# Set run length (keep it to 18 hours for stability if needed, or 48 for full run)
init_hour = int(H_init.date.strftime('%H'))
max_fxx = 48 if init_hour in [0, 6, 12, 18] else 18

print(f"Starting generation for {max_fxx} frames...")

# 2. Loop through frames
for fxx in range(max_fxx + 1):
    try:
        # Re-initialize Herbie for the specific forecast hour
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), 
                   model='hrrr', product='sfc', fxx=fxx, priority=['aws'])

        # Get the projection from the Herbie object directly (Reliable)
        map_crs = H.crs
        
        # Get time strings
        utc_v = H.valid_date.replace(tzinfo=pytz.UTC)
        et_v = utc_v.astimezone(pytz.timezone('US/Eastern'))
        timestamp = f"Valid: {utc_v.strftime('%m/%d %H:%M')}Z | {et_v.strftime('%I:%M %p ET')}"

        # --- PART A: TEMPERATURE MAP ---
        try:
            # Fetch ONLY temperature
            ds_t = H.xarray("TMP:2 m")
            
            # Safety: If it's a list, take the first item
            if isinstance(ds_t, list):
                ds_t = ds_t[0]

            temp_f = (ds_t.t2m - 273.15) * 9/5 + 32

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

        # --- PART B: PRECIP TYPE MAP ---
        try:
            # Fetch ONLY precip types
            ds_p = H.xarray(":(CRAIN|CSNOW|CFRZR|CICEP):")
            
            # Initialize mask using the temperature grid shape (safe bet)
            precip_mask = np.zeros_like(temp_f)

            # Helper function to process a dataset
            def process_precip_ds(d):
                mask = np.zeros_like(d.t2m if 't2m' in d else list(d.data_vars.values())[0])
                if 'crain' in d: mask = np.where(d.crain == 1, 1, mask)
                if 'csnow' in d: mask = np.where(d.csnow == 1, 2, mask)
                if 'cfrzr' in d: mask = np.where(d.cfrzr == 1, 3, mask)
                if 'cicep' in d: mask = np.where(d.cicep == 1, 4, mask)
                return mask

            # Handle List vs Single Object
            if isinstance(ds_p, list):
                for sub_ds in ds_p:
                    # We add the masks together (simple approach)
                    precip_mask = np.maximum(precip_mask, process_precip_ds(sub_ds))
            else:
                precip_mask = process_precip_ds(ds_p)

            fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection': map_crs})
            ax.set_extent(ZOOM_BOUNDS, crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
            
            # Custom Color Map: Transparent, Rain(Gn), Snow(Bl), Ice(Rd), Mix(Or)
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
