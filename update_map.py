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
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

# --- CONFIGURATION ---
EXTENTS = [-80.5, -71.5, 38.5, 43.5] 
MAP_CRS = ccrs.LambertConformal(central_longitude=-76.0, central_latitude=41.0)

# Pre-load County Shapes for background context
try:
    reader = shpreader.Reader(shpreader.natural_earth(resolution='10m', category='cultural', name='admin_2_counties'))
    counties = list(reader.records())
except:
    counties = []

# --- COLOR PALETTES ---

# Wind & Gust (Screenshot: image_ec51d1.png)
WIND_COLORS = [
    '#1e466e', '#2c69b0', '#4292c6', '#6baed6', '#9ecae1', # 0-40
    '#c6dbef', '#e5f5e0', '#a1d99b', '#74c476', '#31a354', # 40-80
    '#fd8d3c', '#f03b20', '#bd0026', '#800026', '#49006a', # 80-120
    '#f768a1', '#fa9fb5', '#fcc5c0', '#fee0d2', '#ffffff'  # 120-170+
]
WIND_LEVELS = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 200]

# Temp & Wind Chill (Screenshot: image_ec4e17.png)
TEMP_COLORS = [
    '#7a0177', '#ae017e', '#dd3497', '#f768a1', '#fbb4b9', # -60 to -20
    '#fee0d2', '#ffffff', '#e0f3f8', '#abd9e9', '#74add1', # -20 to 20
    '#4575b4', '#313695', '#00441b', '#006d2c', '#238b45', # 20 to 60
    '#41ab5d', '#74c476', '#a1d99b', '#e5f5e0', '#ffffcc', # 60 to 80
    '#fed976', '#feb24c', '#fd8d3c', '#f03b20', '#bd0026', # 80 to 110
    '#800026', '#49006a'                                  # 110 to 120+
]
TEMP_LEVELS = list(range(-60, 130, 5))

# Total Precip (No White - Gap filled with Lime/Green)
PRECIP_COLORS = [
    '#ccff99', '#99ff33', '#00cc00', '#006600', '#004d66', # 0.01 to 1.0
    '#3399ff', '#00ffff', '#9999ff', '#9933ff', '#cc33ff', # 1.0 to 2.5
    '#990000', '#cc0000', '#ff3300', '#ff9900', '#cc6600'  # 3.0 to 50.0
]
PRECIP_LEVELS = [0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, 50.0]

# --- HELPERS ---

def get_data_array(ds):
    if isinstance(ds, list): ds = xr.merge(ds, compat='override')
    vars = [v for v in ds.data_vars]
    return ds, ds[vars[0]]

def setup_map(title):
    fig = plt.figure(figsize=(12, 10), facecolor='#020617')
    ax = plt.axes(projection=MAP_CRS)
    ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='black', linewidth=1.5, zorder=3)
    plt.title(title, loc='left', fontweight='bold', fontsize=12, color='white')
    return fig, ax

# --- INITIALIZATION ---
H_init = None
print("Searching for latest HRRR run...")
for offset in range(24):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        test_H = Herbie(search_time.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=1)
        if test_H.inventory() is not None:
            H_init = test_H
            break
    except: continue

if not H_init: 
    print("Could not find data.")
    exit(1)

folders = ["frames_temp", "frames_chill", "frames_precip", "frames_wind", "frames_gust", "frames_snow", "frames_total_precip"]
for folder in folders: os.makedirs(folder, exist_ok=True)

# --- GENERATION LOOP ---
for fxx in range(1, 19):
    try:
        H = Herbie(H_init.date, model='hrrr', fxx=fxx)
        valid_t = H.valid_date.strftime('%m/%d %I:%M %p')

        # 1. TEMPERATURE
        ds_raw = H.xarray("TMP:2 m", verbose=False)
        ds, temp_k = get_data_array(ds_raw)
        data = (temp_k - 273.15) * 9/5 + 32
        fig, ax = setup_map(f"HRRR 2m Temp (°F) | {valid_t}")
        ax.pcolormesh(ds.longitude, ds.latitude, data, transform=ccrs.PlateCarree(), 
                      cmap=ListedColormap(TEMP_COLORS), norm=BoundaryNorm(TEMP_LEVELS, len(TEMP_COLORS)))
        plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()

        # 2. WIND CHILL
        try:
            ds_raw = H.xarray(":WCHILL:surface", verbose=False)
            ds, chill_k = get_data_array(ds_raw)
            data = (chill_k - 273.15) * 9/5 + 32
            fig, ax = setup_map(f"HRRR Wind Chill (°F) | {valid_t}")
            ax.pcolormesh(ds.longitude, ds.latitude, data, transform=ccrs.PlateCarree(), 
                          cmap=ListedColormap(TEMP_COLORS), norm=BoundaryNorm(TEMP_LEVELS, len(TEMP_COLORS)))
            plt.savefig(f"frames_chill/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except: pass

        # 3. WIND SPEED (Separate Section)
        try:
            ds_w = H.xarray(":(UGRD|VGRD):10 m", verbose=False)
            ds, _ = get_data_array(ds_w)
            ws = np.sqrt(ds['u10']**2 + ds['v10']**2) * 2.237 # knots to mph
            fig, ax = setup_map(f"HRRR Wind Speed (mph) | {valid_t}")
            ax.pcolormesh(ds.longitude, ds.latitude, ws, transform=ccrs.PlateCarree(), 
                          cmap=ListedColormap(WIND_COLORS), norm=BoundaryNorm(WIND_LEVELS, len(WIND_COLORS)))
            plt.savefig(f"frames_wind/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except: pass

        # 4. WIND GUSTS (Separate Section)
        try:
            ds_raw = H.xarray(":GUST:surface", verbose=False)
            ds, var = get_data_array(ds_raw)
            gusts = var * 2.237
            fig, ax = setup_map(f"HRRR Wind Gust (mph) | {valid_t}")
            ax.pcolormesh(ds.longitude, ds.latitude, gusts, transform=ccrs.PlateCarree(), 
                          cmap=ListedColormap(WIND_COLORS), norm=BoundaryNorm(WIND_LEVELS, len(WIND_COLORS)))
            plt.savefig(f"frames_gust/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except: pass

        # 5. PRECIP TYPE (Darker Snow Blue)
        try:
            ds_p = H.xarray(":C(RAIN|FREEZ|ICEP|SNOW):surface", verbose=False)
            # 1=Rain(Green), 2=Mix(Orange), 3=Ice(Red), 4=Snow(Darker Blue)
            ptype_data = ds_p.crain*1 + ds_p.cfrzr*2 + ds_p.icep*3 + ds_p.csnow*4
            fig, ax = setup_map(f"HRRR Precip Type | {valid_t}")
            ax.pcolormesh(ds_p.longitude, ds_p.latitude, ptype_data.where(ptype_data > 0), 
                          transform=ccrs.PlateCarree(), 
                          cmap=ListedColormap(['#33ff33', '#ff9900', '#ff0000', '#0044cc']), # Darker Blue for Snow
                          norm=BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5], 4))
            plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except: pass

        # 6. TOTAL PRECIP (No White)
        try:
            ds_raw = H.xarray(":APCP:surface", verbose=False)
            ds, var = get_data_array(ds_raw)
            precip_in = var * 0.03937
            fig, ax = setup_map(f"HRRR Total Precip (in) | {valid_t}")
            # Mask out 0 to make it clear, starting from 0.01 with lime color
            ax.pcolormesh(ds.longitude, ds.latitude, precip_in.where(precip_in >= 0.01), 
                          transform=ccrs.PlateCarree(), cmap=ListedColormap(PRECIP_COLORS), 
                          norm=BoundaryNorm(PRECIP_LEVELS, len(PRECIP_COLORS)))
            plt.savefig(f"frames_total_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except: pass

        # 7. TOTAL SNOWFALL
        try:
            ds_raw = H.xarray(":ASNOW:surface", verbose=False)
            ds, snow_m = get_data_array(ds_raw)
            snow_in = snow_m * 39.37
            fig, ax = setup_map(f"HRRR Total Snowfall (in) | {valid_t}")
            ax.pcolormesh(ds.longitude, ds.latitude, snow_in.where(snow_in >= 0.1), transform=ccrs.PlateCarree(), 
                          cmap='Blues', vmin=0, vmax=12) # Standard blue scale for snow depth
            plt.savefig(f"frames_snow/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except: pass

        print(f"Completed Frame f{fxx:02d}")

    except Exception as e:
        print(f"Fatal error on f{fxx}: {e}")
