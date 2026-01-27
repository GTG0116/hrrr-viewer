import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from herbie import Herbie
from datetime import datetime, timedelta
import pytz
import os
import numpy as np
import xarray as xr
from matplotlib.colors import ListedColormap, BoundaryNorm
import warnings

warnings.filterwarnings("ignore")

# --- CONFIGURATION & STYLE ---
ZOOM_BOUNDS = [-82, -67, 37, 48] 
BG_COLOR = "#0b131e"  # Matches your screenshot's dark blue
TEXT_COLOR = "#ffffff"

# 1. Find Data
H_init = None
for offset in range(24):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        date_str = search_time.strftime('%Y-%m-%d %H:00')
        test_H = Herbie(date_str, model='hrrr', product='sfc', fxx=0, priority=['aws'])
        if test_H.inventory() is not None:
            H_init = test_H
            break
    except: continue

if not H_init: raise Exception("Data unavailable.")

os.makedirs("frames_temp", exist_ok=True)
os.makedirs("frames_precip", exist_ok=True)

# 2. Plotting Loop
for fxx in range(19): # Generating first 18 hours for speed
    try:
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), model='hrrr', product='sfc', fxx=fxx)
        utc_v = H.valid_date.replace(tzinfo=pytz.UTC)
        et_v = utc_v.astimezone(pytz.timezone('US/Eastern'))
        timestamp = f"{et_v.strftime('%I:%M %p ET')} | {et_v.strftime('%a %b %d')}"

        # Get Temp and Projection
        ds_t = H.xarray("TMP:2 m")
        if isinstance(ds_t, list): ds_t = ds_t[0]
        temp_f = (ds_t.t2m - 273.15) * 9/5 + 32
        
        try: map_crs = ds_t.herbie.crs
        except: map_crs = ccrs.LambertConformal(central_longitude=-97.5, central_latitude=38.5)

        # --- MAP 1: TEMPERATURE (Sleek Gradient) ---
        fig = plt.figure(figsize=(10, 8), facecolor=BG_COLOR)
        ax = plt.axes(projection=map_crs, facecolor=BG_COLOR)
        ax.set_extent(ZOOM_BOUNDS, crs=ccrs.PlateCarree())
        
        # Style the map
        ax.add_feature(cfeature.STATES, edgecolor='#444444', linewidth=0.8)
        ax.add_feature(cfeature.COASTLINE, edgecolor='#444444', linewidth=0.8)
        
        im = ax.pcolormesh(ds_t.longitude, ds_t.latitude, temp_f, transform=ccrs.PlateCarree(), 
                           cmap='magma', vmin=0, vmax=100, alpha=0.8)
        
        # Text Overlays
        plt.text(0.02, 0.94, "EPHRATA WEATHER", transform=ax.transAxes, color=TEXT_COLOR, 
                 fontsize=14, fontweight='bold', family='sans-serif')
        plt.text(0.02, 0.90, f"HRRR 2m Temperature • {timestamp}", transform=ax.transAxes, 
                 color='#aaaaaa', fontsize=10)
        
        plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor=BG_COLOR)
        plt.close()

        # --- MAP 2: PRECIP TYPE (Intensity Shaded) ---
        # Fetching simulated reflectivity for intensity
        ds_p = H.xarray(":(CRAIN|CSNOW|CFRZR|CICEP|REFC):")
        if isinstance(ds_p, list): ds_p = xr.merge(ds_p)
        
        refc = ds_p.refc.values if 'refc' in ds_p else np.zeros_like(temp_f)
        
        # Logic for shaded precip types
        # 1-10: Rain (Light to Heavy), 11-20: Snow, 21-30: Ice
        shaded_mask = np.zeros_like(temp_f)
        if 'crain' in ds_p: shaded_mask = np.where(ds_p.crain == 1, np.clip(refc/5, 1, 9), shaded_mask)
        if 'csnow' in ds_p: shaded_mask = np.where(ds_p.csnow == 1, np.clip(refc/5 + 10, 11, 19), shaded_mask)
        
        # Custom "Vibe" Colormap
        # Shades of Green (Rain), Blues (Snow), Pinks (Ice)
        colors = ['#0b131e'] # Transparent/BG
        colors.extend(['#1a4314', '#2e7d32', '#4caf50', '#81c784']) # Green rain
        colors.extend(['#0d47a1', '#1976d2', '#42a5f5', '#90caf9']) # Blue snow
        p_cmap = ListedColormap(colors)

        fig = plt.figure(figsize=(10, 8), facecolor=BG_COLOR)
        ax = plt.axes(projection=map_crs, facecolor=BG_COLOR)
        ax.set_extent(ZOOM_BOUNDS, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.STATES, edgecolor='#555555', linewidth=1)
        
        ax.pcolormesh(ds_t.longitude, ds_t.latitude, np.ma.masked_where(shaded_mask == 0, shaded_mask),
                      transform=ccrs.PlateCarree(), cmap=p_cmap)

        # Dashboard Title
        plt.text(0.02, 0.94, "EPHRATA WEATHER", transform=ax.transAxes, color=TEXT_COLOR, 
                 fontsize=14, fontweight='bold')
        plt.text(0.02, 0.90, f"High-Res Forecast Radar • {timestamp}", transform=ax.transAxes, color='#aaaaaa')
        
        # Legend (Manual)
        plt.text(0.80, 0.05, "RAIN", color='#4caf50', transform=ax.transAxes, fontweight='bold')
        plt.text(0.88, 0.05, "SNOW", color='#42a5f5', transform=ax.transAxes, fontweight='bold')

        plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor=BG_COLOR)
        plt.close()
        print(f"Rendered f{fxx}")

    except Exception as e: print(f"Error f{fxx}: {e}")
