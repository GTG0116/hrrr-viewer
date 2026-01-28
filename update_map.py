import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from herbie import Herbie
from datetime import datetime, timedelta
import os
import numpy as np
import xarray as xr
from matplotlib.colors import ListedColormap, BoundaryNorm
import warnings

warnings.filterwarnings("ignore")

# --- CONFIGURATION ---
EXTENTS = [-80.5, -71.5, 38.5, 43.5] 
MAP_CRS = ccrs.LambertConformal(central_longitude=-76.0, central_latitude=41.0)

# --- CORRECTED COLOR PALETTES ---

# Wind & Gust Palette (Matches image_ec51d1.png)
# 18 Levels require 17 Colors
WIND_COLORS = [
    '#1e466e', '#2c69b0', '#4292c6', '#6baed6', '#9ecae1', 
    '#c6dbef', '#e5f5e0', '#a1d99b', '#74c476', '#31a354', 
    '#fd8d3c', '#f03b20', '#bd0026', '#800026', '#49006a', 
    '#f768a1', '#fa9fb5'
]
WIND_LEVELS = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170]

# Temperature & Wind Chill (Matches image_ec4e17.png)
# 37 Levels require 36 Colors
TEMP_COLORS = [
    '#7a0177', '#ae017e', '#dd3497', '#f768a1', '#fbb4b9', '#fee0d2', '#ffffff', 
    '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695', '#00441b', '#006d2c', 
    '#238b45', '#41ab5d', '#74c476', '#a1d99b', '#e5f5e0', '#ffffcc', '#fed976', 
    '#feb24c', '#fd8d3c', '#f03b20', '#bd0026', '#800026', '#49006a', '#4d004b',
    '#7b3294', '#c2a5cf', '#f7f7f7', '#a6dba0', '#008837', '#fee08b', '#d73027', '#67001f'
]
TEMP_LEVELS = list(range(-60, 121, 5))

# Total Precip Palette (Matches Screenshot 2026-01-27 215251.png)
# Removed white; gap between 0 and 0.1 is filled with lime (#ccff99)
PRECIP_COLORS = [
    '#ccff99', '#99ff33', '#00cc00', '#006600', '#004d66', 
    '#3399ff', '#00ffff', '#9999ff', '#9933ff', '#cc33ff', 
    '#990000', '#cc0000', '#ff3300', '#ff9900', '#cc6600'
]
PRECIP_LEVELS = [0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, 20.0]

# --- HELPERS ---

def setup_map(title):
    fig = plt.figure(figsize=(12, 10), facecolor='#020617')
    ax = plt.axes(projection=MAP_CRS)
    ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='black', linewidth=1.5, zorder=3)
    plt.title(title, loc='left', fontweight='bold', fontsize=12, color='white')
    return fig, ax

# --- INITIALIZATION (Fixes H_init NameError) ---
H_init = None
for offset in range(12):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        test_H = Herbie(search_time.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=1)
        if test_H.inventory() is not None:
            H_init = test_H
            break
    except: continue

if not H_init: 
    print("Error: Could not find recent HRRR data.")
    exit(1)

folders = ["frames_temp", "frames_chill", "frames_precip", "frames_wind", "frames_gust", "frames_snow", "frames_total_precip"]
for f in folders: os.makedirs(f, exist_ok=True)

# --- GENERATION LOOP ---
for fxx in range(1, 19):
    try:
        H = Herbie(H_init.date, model='hrrr', fxx=fxx)
        valid_t = H.valid_date.strftime('%m/%d %I:%M %p')

        # 1. TEMPERATURE & WIND CHILL
        for param, folder, name in [("TMP:2 m", "frames_temp", "Temperature"), (":WCHILL:", "frames_chill", "Wind Chill")]:
            ds = H.xarray(param, verbose=False)
            data = (ds[list(ds.data_vars)[0]] - 273.15) * 9/5 + 32
            fig, ax = setup_map(f"HRRR {name} (°F) | {valid_t}")
            ax.pcolormesh(ds.longitude, ds.latitude, data, transform=ccrs.PlateCarree(), 
                          cmap=ListedColormap(TEMP_COLORS), norm=BoundaryNorm(TEMP_LEVELS, len(TEMP_COLORS)))
            plt.savefig(f"{folder}/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()

        # 2. WIND SPEED
        ds_w = H.xarray(":(UGRD|VGRD):10 m", verbose=False)
        ws = np.sqrt(ds_w['u10']**2 + ds_w['v10']**2) * 2.237
        fig, ax = setup_map(f"HRRR Wind Speed (mph) | {valid_t}")
        ax.pcolormesh(ds_w.longitude, ds_w.latitude, ws, transform=ccrs.PlateCarree(), 
                      cmap=ListedColormap(WIND_COLORS), norm=BoundaryNorm(WIND_LEVELS, len(WIND_COLORS)))
        plt.savefig(f"frames_wind/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()

        # 3. WIND GUSTS
        ds_g = H.xarray(":GUST:surface", verbose=False)
        gusts = ds_g['gust'] * 2.237
        fig, ax = setup_map(f"HRRR Wind Gust (mph) | {valid_t}")
        ax.pcolormesh(ds_g.longitude, ds_g.latitude, gusts, transform=ccrs.PlateCarree(), 
                      cmap=ListedColormap(WIND_COLORS), norm=BoundaryNorm(WIND_LEVELS, len(WIND_COLORS)))
        plt.savefig(f"frames_gust/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()

        # 4. PRECIP TYPE (Darker Snow Blue)
        ds_p = H.xarray(":C(RAIN|FREEZ|ICEP|SNOW):surface", verbose=False)
        ptype = ds_p.crain*1 + ds_p.cfrzr*2 + ds_p.icep*3 + ds_p.csnow*4
        fig, ax = setup_map(f"HRRR Precip Type | {valid_t}")
        # Snow category updated to a deep midnight blue (#001f3f)
        ax.pcolormesh(ds_p.longitude, ds_p.latitude, ptype.where(ptype > 0), transform=ccrs.PlateCarree(), 
                      cmap=ListedColormap(['#33ff33', '#ff9900', '#ff0000', '#001f3f']), 
                      norm=BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5], 4))
        plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()

        # 5. TOTAL PRECIP (No White)
        ds_tp = H.xarray(":APCP:surface", verbose=False)
        precip = ds_tp['tp'] * 0.03937
        fig, ax = setup_map(f"HRRR Total Precip (in) | {valid_t}")
        ax.pcolormesh(ds_tp.longitude, ds_tp.latitude, precip.where(precip >= 0.01), transform=ccrs.PlateCarree(), 
                      cmap=ListedColormap(PRECIP_COLORS), norm=BoundaryNorm(PRECIP_LEVELS, len(PRECIP_COLORS)))
        plt.savefig(f"frames_total_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()

        print(f"Completed Frame f{fxx:02d}")
    except Exception as e:
        print(f"Error on frame f{fxx:02d}: {e}")
