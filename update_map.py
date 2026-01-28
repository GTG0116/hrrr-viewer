import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from herbie import Herbie
from datetime import datetime, timedelta
import pytz
import os
import numpy as np
import xarray as xr
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.patches as mpatches
import warnings

# Suppress annoying warnings
warnings.filterwarnings("ignore")

# --- CONFIGURATION ---
EXTENTS = [-80.5, -71.5, 38.5, 43.5] 
MAP_CRS = ccrs.LambertConformal(central_longitude=-76.0, central_latitude=41.0)

# Pre-load County Shapes
try:
    reader = shpreader.Reader(shpreader.natural_earth(resolution='10m', category='cultural', name='admin_2_counties'))
    counties = list(reader.records())
except:
    counties = []

# --- COLOR PALETTES (FROM SCREENSHOTS) ---

# Wind & Gust Palette (image_ec51d1.png)
WIND_COLORS = [
    '#1e466e', '#2c69b0', '#4292c6', '#6baed6', '#9ecae1', # 0-40
    '#c6dbef', '#e5f5e0', '#a1d99b', '#74c476', '#31a354', # 40-80
    '#fd8d3c', '#f03b20', '#bd0026', '#800026', '#49006a', # 80-120
    '#f768a1', '#fa9fb5', '#fcc5c0', '#fee0d2', '#ffffff'  # 120-170+
]
WIND_LEVELS = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 200]

# Temperature & Wind Chill Palette (image_ec4e17.png)
TEMP_COLORS = [
    '#7a0177', '#ae017e', '#dd3497', '#f768a1', '#fbb4b9', # -60 to -20
    '#fee0d2', '#ffffff', '#e0f3f8', '#abd9e9', '#74add1', # -20 to 20
    '#4575b4', '#313695', '#00441b', '#006d2c', '#238b45', # 20 to 60
    '#41ab5d', '#74c476', '#a1d99b', '#e5f5e0', '#ffffcc', # 60 to 80
    '#fed976', '#feb24c', '#fd8d3c', '#f03b20', '#bd0026', # 80 to 110
    '#800026', '#49006a'                                  # 110 to 120+
]
TEMP_LEVELS = list(range(-60, 130, 5)) # Detailed 5-degree increments

# Total Precip Palette (No White - Gap filled with light green/blue)
PRECIP_COLORS = [
    '#ffffff00', # Transparent (0)
    '#ccff99', '#99ff33', '#00cc00', '#006600', '#004d66', # 0.01 to 1.0
    '#3399ff', '#00ffff', '#9999ff', '#9933ff', '#cc33ff', # 1.0 to 2.5
    '#990000', '#cc0000', '#ff3300', '#ff9900', '#cc6600'  # 3.0 to 50.0
]
PRECIP_LEVELS = [0, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, 50.0]

# --- HELPERS ---

def get_data_array(ds):
    if isinstance(ds, list): ds = xr.merge(ds, compat='override')
    vars = [v for v in ds.data_vars]
    if not vars: raise ValueError("No data variables found")
    return ds, ds[vars[0]]

def setup_map(title):
    fig = plt.figure(figsize=(12, 10), facecolor='#020617')
    ax = plt.axes(projection=MAP_CRS)
    ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='black', linewidth=1.2, zorder=2)
    ax.add_feature(cfeature.ShapelyFeature([c.geometry for c in counties], ccrs.PlateCarree()), 
                   facecolor='none', edgecolor='gray', linewidth=0.4, alpha=0.6, zorder=1)
    plt.title(title, loc='left', fontweight='bold', fontsize=12, color='white')
    return fig, ax

# --- DATA SEARCH & MAIN LOOP ---
# (Search logic remains the same as previous version)
# ... [Herbie Search Logic] ...

for fxx in range(1, 19):
    H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=fxx)
    t_str = f"Valid: {H.valid_date.strftime('%m/%d %I:%M %p')}"

    # 1. TEMPERATURE & WIND CHILL
    for param, folder, name in [("TMP:2 m", "frames_temp", "Temperature"), (":WCHILL:", "frames_chill", "Wind Chill")]:
        try:
            ds_raw = H.xarray(param, verbose=False)
            ds, var = get_data_array(ds_raw)
            data_f = (var - 273.15) * 9/5 + 32
            fig, ax = setup_map(f"HRRR 2m {name} (°F) | {t_str}")
            im = ax.pcolormesh(ds.longitude, ds.latitude, data_f, transform=ccrs.PlateCarree(), 
                               cmap=ListedColormap(TEMP_COLORS), norm=BoundaryNorm(TEMP_LEVELS, len(TEMP_COLORS)))
            plt.savefig(f"{folder}/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except: continue

    # 2. WIND & GUSTS (Using Wind Palette)
    for param, folder, name in [(":(UGRD|VGRD):10 m", "frames_wind", "Wind Speed"), (":GUST:surface", "frames_gust", "Wind Gust")]:
        try:
            ds_raw = H.xarray(param, verbose=False)
            ds, var = get_data_array(ds_raw)
            # Calculate speed if U/V, otherwise use direct gust var
            speed = (np.sqrt(ds['u10']**2 + ds['v10']**2) * 2.237) if 'u10' in ds else (var * 2.237)
            fig, ax = setup_map(f"HRRR {name} (mph) | {t_str}")
            im = ax.pcolormesh(ds.longitude, ds.latitude, speed, transform=ccrs.PlateCarree(), 
                               cmap=ListedColormap(WIND_COLORS), norm=BoundaryNorm(WIND_LEVELS, len(WIND_COLORS)))
            plt.savefig(f"{folder}/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except: continue

    # 3. TOTAL PRECIP (No White Scale)
    try:
        ds_raw = H.xarray(":APCP:surface", verbose=False)
        ds, var = get_data_array(ds_raw)
        precip_in = var * 0.03937
        fig, ax = setup_map(f"HRRR Total Precip (in) | {t_str}")
        im = ax.pcolormesh(ds.longitude, ds.latitude, precip_in, transform=ccrs.PlateCarree(), 
                           cmap=ListedColormap(PRECIP_COLORS), norm=BoundaryNorm(PRECIP_LEVELS, len(PRECIP_COLORS)))
        plt.savefig(f"frames_total_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()
    except: continue
