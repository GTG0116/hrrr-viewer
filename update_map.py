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
SAVE_DIR = os.path.join(os.getcwd(), "herbie_data") # Force local download

# Pre-load County Shapes
try:
    reader = shpreader.Reader(shpreader.natural_earth(resolution='10m', category='cultural', name='admin_2_counties'))
    counties = list(reader.records())
except:
    counties = []

# --- COLOR PALETTES ---

# Wind & Gust (Blue -> Pink)
WIND_COLORS = [
    '#1e466e', '#2c69b0', '#4292c6', '#6baed6', '#9ecae1', 
    '#c6dbef', '#e5f5e0', '#a1d99b', '#74c476', '#31a354', 
    '#fd8d3c', '#f03b20', '#bd0026', '#800026', '#49006a', 
    '#f768a1', '#fa9fb5'
]
WIND_LEVELS = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170]

# Temp & Chill (-60 to 120)
TEMP_COLORS = [
    '#7a0177', '#ae017e', '#dd3497', '#f768a1', '#fbb4b9', '#fee0d2', '#ffffff', 
    '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695', '#00441b', '#006d2c', 
    '#238b45', '#41ab5d', '#74c476', '#a1d99b', '#e5f5e0', '#ffffcc', '#fed976', 
    '#feb24c', '#fd8d3c', '#f03b20', '#bd0026', '#800026', '#49006a', '#4d004b',
    '#7b3294', '#c2a5cf', '#f7f7f7', '#a6dba0', '#008837', '#fee08b', '#d73027', '#67001f'
]
TEMP_LEVELS = list(range(-60, 121, 5))

# Total Precip (Lime start, no white)
PRECIP_COLORS = [
    '#ccff99', '#99ff33', '#00cc00', '#006600', '#004d66', 
    '#3399ff', '#00ffff', '#9999ff', '#9933ff', '#cc33ff', 
    '#990000', '#cc0000', '#ff3300', '#ff9900', '#cc6600'
]
PRECIP_LEVELS = [0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, 20.0]

# Snowfall (Blues)
SNOW_COLORS = [
    '#d1e5f0', '#92c5de', '#4393c3', '#2166ac', '#053061', 
    '#ffffbf', '#fed976', '#feb24c', '#fd8d3c', '#f03b20', 
    '#bd0026', '#800026', '#49000a', '#d8daeb', '#b2abd2'
]
SNOW_LEVELS = [0.1, 1, 2, 3, 4, 6, 8, 12, 18, 24, 30, 36, 48, 60, 72, 100]

# --- HELPERS ---

def setup_map(title):
    fig = plt.figure(figsize=(12, 10), facecolor='#020617')
    ax = plt.axes(projection=MAP_CRS)
    ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='black', linewidth=1.5, zorder=3)
    plt.title(title, loc='left', fontweight='bold', fontsize=12, color='white')
    return fig, ax

def draw_county_labels(ax, ds, data_var, fmt="{:.0f}", check_val=None):
    """Draws values at county centroids."""
    view_w, view_e, view_s, view_n = EXTENTS
    for county in counties:
        if not county.geometry.centroid.within(county.geometry): continue 
        lon, lat = county.geometry.centroid.x, county.geometry.centroid.y
        if lon < view_w or lon > view_e or lat < view_s or lat > view_n: continue

        try:
            val = float(data_var.sel(latitude=lat, longitude=lon, method='nearest').values)
            if check_val is not None and val < check_val: continue
            ax.text(lon, lat, fmt.format(val), transform=ccrs.PlateCarree(),
                    ha='center', va='center', fontsize=7, fontweight='bold',
                    color='black', clip_on=True, zorder=10)
        except: continue

# --- INITIALIZATION ---
H_init = None
for offset in range(12):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        # Add save_dir here to fix the file not found error
        test_H = Herbie(search_time.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=1, save_dir=SAVE_DIR)
        if test_H.inventory() is not None:
            H_init = test_H
            break
    except: continue

if not H_init: 
    print("Error: Could not find HRRR data.")
    exit(1)

folders = ["frames_temp", "frames_chill", "frames_precip", "frames_wind", "frames_gust", "frames_snow", "frames_total_precip"]
for f in folders: os.makedirs(f, exist_ok=True)

# --- GENERATION LOOP ---
for fxx in range(1, 19):
    try:
        # Add save_dir to the loop Herbie object as well
        H = Herbie(H_init.date, model='hrrr', fxx=fxx, save_dir=SAVE_DIR)
        valid_t = H.valid_date.strftime('%m/%d %I:%M %p')

        # 1. TEMPERATURE & WIND CHILL
        for param, folder, name, check in [("TMP:2 m", "frames_temp", "Temperature", None), (":WCHILL:", "frames_chill", "Wind Chill", None)]:
            try:
                ds = H.xarray(param, verbose=False)
                data = (ds[list(ds.data_vars)[0]] - 273.15) * 9/5 + 32
                fig, ax = setup_map(f"HRRR {name} (°F) | {valid_t}")
                ax.pcolormesh(ds.longitude, ds.latitude, data, transform=ccrs.PlateCarree(), 
                            cmap=ListedColormap(TEMP_COLORS), norm=BoundaryNorm(TEMP_LEVELS, len(TEMP_COLORS)))
                draw_county_labels(ax, ds, data, fmt="{:.0f}")
                plt.savefig(f"{folder}/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
                plt.close()
            except: pass

        # 2. WIND SPEED (Separate)
        try:
            ds_w = H.xarray(":(UGRD|VGRD):10 m", verbose=False)
            ws = np.sqrt(ds_w['u10']**2 + ds_w['v10']**2) * 2.237
            fig, ax = setup_map(f"HRRR Wind Speed (mph) | {valid_t}")
            ax.pcolormesh(ds_w.longitude, ds_w.latitude, ws, transform=ccrs.PlateCarree(), 
                        cmap=ListedColormap(WIND_COLORS), norm=BoundaryNorm(WIND_LEVELS, len(WIND_COLORS)))
            draw_county_labels(ax, ds_w, ws, fmt="{:.0f}", check_val=10)
            plt.savefig(f"frames_wind/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except: pass

        # 3. WIND GUSTS (Separate)
        try:
            ds_g = H.xarray(":GUST:surface", verbose=False)
            gusts = ds_g['gust'] * 2.237
            fig, ax = setup_map(f"HRRR Wind Gust (mph) | {valid_t}")
            ax.pcolormesh(ds_g.longitude, ds_g.latitude, gusts, transform=ccrs.PlateCarree(), 
                        cmap=ListedColormap(WIND_COLORS), norm=BoundaryNorm(WIND_LEVELS, len(WIND_COLORS)))
            draw_county_labels(ax, ds_g, gusts, fmt="{:.0f}", check_val=20)
            plt.savefig(f"frames_gust/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except: pass

        # 4. PRECIP TYPE (Midnight Blue Snow)
        try:
            ds_p = H.xarray(":C(RAIN|FREEZ|ICEP|SNOW):surface", verbose=False)
            ptype = ds_p.crain*1 + ds_p.cfrzr*2 + ds_p.icep*3 + ds_p.csnow*4
            fig, ax = setup_map(f"HRRR Precip Type | {valid_t}")
            ax.pcolormesh(ds_p.longitude, ds_p.latitude, ptype.where(ptype > 0), transform=ccrs.PlateCarree(), 
                        cmap=ListedColormap(['#33ff33', '#ff9900', '#ff0000', '#001f3f']), 
                        norm=BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5], 4))
            plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except: pass

        # 5. TOTAL PRECIP (No White)
        try:
            ds_tp = H.xarray(":APCP:surface", verbose=False)
            precip = ds_tp['tp'] * 0.03937
            fig, ax = setup_map(f"HRRR Total Precip (in) | {valid_t}")
            ax.pcolormesh(ds_tp.longitude, ds_tp.latitude, precip.where(precip >= 0.01), transform=ccrs.PlateCarree(), 
                        cmap=ListedColormap(PRECIP_COLORS), norm=BoundaryNorm(PRECIP_LEVELS, len(PRECIP_COLORS)))
            draw_county_labels(ax, ds_tp, precip, fmt="{:.2f}", check_val=0.05)
            plt.savefig(f"frames_total_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except: pass

        # 6. TOTAL SNOWFALL
        try:
            ds_s = H.xarray(":ASNOW:surface", verbose=False)
            snow = ds_s['asnow'] * 39.37
            fig, ax = setup_map(f"HRRR Total Snow (in) | {valid_t}")
            ax.pcolormesh(ds_s.longitude, ds_s.latitude, snow.where(snow >= 0.1), transform=ccrs.PlateCarree(), 
                        cmap=ListedColormap(SNOW_COLORS), norm=BoundaryNorm(SNOW_LEVELS, len(SNOW_COLORS)))
            draw_county_labels(ax, ds_s, snow, fmt="{:.1f}", check_val=0.5)
            plt.savefig(f"frames_snow/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except: pass

        print(f"Completed Frame f{fxx:02d}")
    except Exception as e:
        print(f"Error on frame f{fxx:02d}: {e}")
