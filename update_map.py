import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from herbie import Herbie
from datetime import datetime, timedelta
import pytz
import os
import numpy as np

# --- CONFIGURATION ---
# Change these to zoom in ([West, East, South, North])
# Example for Northeast US: [-82, -67, 37, 48]
# Set to None for full Lower 48
ZOOM_BOUNDS = [-82, -67, 37, 48] 

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
os.makedirs("frames", exist_ok=True)

# 2. Loop through frames
for fxx in range(max_fxx + 1):
    try:
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), 
                   model='hrrr', product='sfc', fxx=fxx, priority=['aws'])
        
        # Pull Temperature AND Precip Types
        ds = H.xarray("(:TMP:2 m|:CSNOW:|:CICEP:|:CFRZR:|:CRAIN:)")
        temp_f = (ds.t2m - 273.15) * 9/5 + 32

        # Create a single Precip Type mask
        # Priority: Ice > Freezing Rain > Snow > Rain
        precip_mask = np.zeros_like(temp_f)
        if hasattr(ds, 'crain'): precip_mask = np.where(ds.crain == 1, 1, precip_mask)
        if hasattr(ds, 'csnow'): precip_mask = np.where(ds.csnow == 1, 2, precip_mask)
        if hasattr(ds, 'cfrzr'): precip_mask = np.where(ds.cfrzr == 1, 3, precip_mask)
        if hasattr(ds, 'cicep'): precip_mask = np.where(ds.cicep == 1, 4, precip_mask)

        # Plotting
        fig = plt.figure(figsize=(12, 8), facecolor='white')
        ax = plt.axes(projection=ds.herbie.crs)
        
        if ZOOM_BOUNDS:
            ax.set_extent(ZOOM_BOUNDS, crs=ccrs.PlateCarree())

        ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='black', linewidth=0.7)
        ax.add_feature(cfeature.LAKES.with_scale('50m'), edgecolor='blue', linewidth=0.5, alpha=0.3)

        # Background: Temperature
        temp_plot = ax.pcolormesh(ds.longitude, ds.latitude, temp_f, 
                                 transform=ccrs.PlateCarree(), cmap='RdYlBu_r', 
                                 vmin=0, vmax=100, alpha=0.8)
        
        # Overlay: Precip Type (using a custom discrete map)
        # 1=Rain(Green), 2=Snow(Blue), 3=FRZR(Red), 4=Ice(Orange)
        from matplotlib.colors import ListedColormap
        precip_cmap = ListedColormap(['none', '#2ecc71', '#3498db', '#e74c3c', '#e67e22'])
        ax.pcolormesh(ds.longitude, ds.latitude, np.ma.masked_where(precip_mask == 0, precip_mask),
                     transform=ccrs.PlateCarree(), cmap=precip_cmap, vmin=0, vmax=4)

        # Labels & Time
        utc_v = H.valid_date.replace(tzinfo=pytz.UTC)
        et_v = utc_v.astimezone(pytz.timezone('US/Eastern'))
        
        plt.title(f"HRRR Temp + Precip Type | Init: {H.date.strftime('%H:%M')}Z\n"
                  f"Valid: {utc_v.strftime('%m/%d %H:%M')}Z ({et_v.strftime('%I:%M %p ET')})", 
                  loc='left', fontweight='bold')

        plt.savefig(f"frames/frame_{fxx:02d}.png", dpi=100, bbox_inches='tight')
        plt.close()
        print(f"Done f{fxx:02d}")

    except Exception as e:
        print(f"Error f{fxx}: {e}")
