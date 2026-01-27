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
ZOOM_BOUNDS = [-82, -67, 37, 48] # [West, East, South, North]

# 1. Find latest data
H_init = None
for offset in range(4):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        date_str = search_time.strftime('%Y-%m-%d %H:00')
        test_H = Herbie(date_str, model='hrrr', product='sfc', fxx=0, priority=['aws'])
        if test_H.inventory() is not None:
            H_init = test_H
            break
    except:
        continue

if not H_init: raise Exception("No HRRR data found.")

init_hour = int(H_init.date.strftime('%H'))
max_fxx = 48 if init_hour in [0, 6, 12, 18] else 18

os.makedirs("frames_temp", exist_ok=True)
os.makedirs("frames_precip", exist_ok=True)

# 2. Loop through frames
for fxx in range(max_fxx + 1):
    try:
        print(f"Processing Hour f{fxx:02d}...")
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), 
                   model='hrrr', product='sfc', fxx=fxx, priority=['aws'])
        
        # Robust Data Loading: Pull Temperature separately from Precip
        ds_temp = H.xarray("TMP:2 m")
        # If it returns a list, take the first item
        if isinstance(ds_temp, list): ds_temp = ds_temp[0]
        
        temp_f = (ds_temp.t2m - 273.15) * 9/5 + 32
        
        # Setup Times
        utc_v = H.valid_date.replace(tzinfo=pytz.UTC)
        et_v = utc_v.astimezone(pytz.timezone('US/Eastern'))
        timestamp = f"Valid: {utc_v.strftime('%m/%d %H:%M')}Z ({et_v.strftime('%I:%M %p ET')})"

        # --- MAP 1: TEMPERATURE ---
        fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection': ds_temp.herbie.crs})
        if ZOOM_BOUNDS: ax.set_extent(ZOOM_BOUNDS, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        im = ax.pcolormesh(ds_temp.longitude, ds_temp.latitude, temp_f, transform=ccrs.PlateCarree(), 
                          cmap='RdYlBu_r', vmin=0, vmax=100)
        plt.colorbar(im, label="Temperature (°F)", orientation='horizontal', pad=0.05)
        plt.title(f"HRRR 2m Temperature\n{timestamp}", loc='left', fontweight='bold')
        plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=100, bbox_inches='tight')
        plt.close()

        # --- MAP 2: PRECIP TYPE ---
        # Pull categorical precip types (CRAIN, CSNOW, CFRZR, CICEP)
        precip_mask = np.zeros_like(temp_f)
        try:
            # We try to get the categorical dataset
            ds_p = H.xarray(":(CRAIN|CSNOW|CFRZR|CICEP):")
            # Handle list return
            if isinstance(ds_p, list):
                # Merge the list into one dataset or process individually
                for sub_ds in ds_p:
                    if 'crain' in sub_ds: precip_mask = np.where(sub_ds.crain == 1, 1, precip_mask)
                    if 'csnow' in sub_ds: precip_mask = np.where(sub_ds.csnow == 1, 2, precip_mask)
                    if 'cfrzr' in sub_ds: precip_mask = np.where(sub_ds.cfrzr == 1, 3, precip_mask)
                    if 'cicep' in sub_ds: precip_mask = np.where(sub_ds.cicep == 1, 4, precip_mask)
        except Exception as pe:
            print(f"No precip found for f{fxx}: {pe}")

        fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection': ds_temp.herbie.crs})
        if ZOOM_BOUNDS: ax.set_extent(ZOOM_BOUNDS, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        p_cmap = ListedColormap(['none', '#2ecc71', '#3498db', '#e74c3c', '#e67e22'])
        ax.pcolormesh(ds_temp.longitude, ds_temp.latitude, np.ma.masked_where(precip_mask == 0, precip_mask),
                     transform=ccrs.PlateCarree(), cmap=p_cmap, vmin=0, vmax=4)
        
        plt.title(f"HRRR Precip Type (Green:Rain, Blue:Snow, Red:Ice)\n{timestamp}", loc='left', fontweight='bold')
        plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=100, bbox_inches='tight')
        plt.close()

    except Exception as e:
        print(f"Critical error on f{fxx}: {e}")
