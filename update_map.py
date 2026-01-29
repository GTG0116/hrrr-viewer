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

warnings.filterwarnings("ignore")

# --- CONFIGURATION ---
EXTENTS = [-80.5, -71.5, 38.5, 43.5] 
MAP_CRS = ccrs.LambertConformal(central_longitude=-76.0, central_latitude=41.0)
SAVE_DIR = os.path.join(os.getcwd(), "herbie_data")

# --- COLOR PALETTES ---
WIND_COLORS = ['#1e466e', '#2c69b0', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef', '#e5f5e0', '#a1d99b', '#74c476', '#31a354', '#fd8d3c', '#f03b20', '#bd0026', '#800026', '#49006a', '#f768a1', '#fa9fb5']
WIND_LEVELS = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170]

TEMP_COLORS = ['#7a0177', '#ae017e', '#dd3497', '#f768a1', '#fbb4b9', '#fee0d2', '#ffffff', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695', '#00441b', '#006d2c', '#238b45', '#41ab5d', '#74c476', '#a1d99b', '#e5f5e0', '#ffffcc', '#fed976', '#feb24c', '#fd8d3c', '#f03b20', '#bd0026', '#800026', '#49006a', '#4d004b', '#7b3294', '#c2a5cf', '#f7f7f7', '#a6dba0', '#008837', '#fee08b', '#d73027', '#67001f']
TEMP_LEVELS = list(range(-60, 121, 5))

PRECIP_COLORS = ['#ccff99', '#99ff33', '#00cc00', '#006600', '#004d66', '#3399ff', '#00ffff', '#9999ff', '#9933ff', '#cc33ff', '#990000', '#cc0000', '#ff3300', '#ff9900', '#cc6600']
PRECIP_LEVELS = [0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, 20.0]

SNOW_COLORS = ['#d1e5f0', '#92c5de', '#4393c3', '#2166ac', '#053061', '#ffffbf', '#fed976', '#feb24c', '#fd8d3c', '#f03b20', '#bd0026', '#800026', '#49000a', '#d8daeb', '#b2abd2']
SNOW_LEVELS = [0.1, 1, 2, 3, 4, 6, 8, 12, 18, 24, 30, 36, 48, 60, 72, 100]

# Pre-load County Shapes
try:
    reader = shpreader.Reader(shpreader.natural_earth(resolution='10m', category='cultural', name='admin_2_counties'))
    counties = list(reader.records())
except: counties = []

# --- HELPERS ---

def safe_get_ds(H, search_str):
    """Retrieves dataset and forces merge if Herbie returns a list (Fixes 'no attribute' errors)."""
    ds = H.xarray(search_str, verbose=False)
    if isinstance(ds, list):
        return xr.merge(ds, compat='override')
    return ds

def draw_county_labels(ax, ds, data_var, fmt="{:.0f}", check_val=None):
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

def setup_plot(H, fxx, datatype):
    """Sets up the figure with date including Year and Ephrata branding."""
    fig = plt.figure(figsize=(12, 10), facecolor='#020617')
    ax = plt.axes(projection=MAP_CRS)
    ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='black', linewidth=1.5, zorder=3)
    
    # Time Conversion
    init_dt = H.date.replace(tzinfo=pytz.UTC)
    valid_dt = H.valid_date.replace(tzinfo=pytz.UTC)
    valid_et = valid_dt.astimezone(pytz.timezone('US/Eastern'))
    
    # Top Left: HRRR: Data Type
    plt.text(0, 1.02, f"HRRR: {datatype}", transform=ax.transAxes, 
             fontsize=14, color='white', fontweight='bold', ha='left')
    
    # Top Right: Run Info (With Year)
    run_str = f"Run: {init_dt.strftime('%m/%d/%Y %H')}Z | F{fxx:02d}"
    valid_str = f"Valid: {valid_dt.strftime('%m/%d/%Y %H')}Z / {valid_et.strftime('%m/%d %I:%M %p')} ET"
    plt.text(1, 1.02, f"{run_str}\n{valid_str}", transform=ax.transAxes, 
             fontsize=9, color='white', ha='right', weight='bold')

    # Bottom Branding
    plt.text(0.5, -0.08, "Ephrata Weather", transform=ax.transAxes,
             fontsize=16, color='white', fontweight='bold', ha='center', va='top')

    return fig, ax

def save_plot(fig, folder, fxx):
    os.makedirs(folder, exist_ok=True)
    plt.savefig(f"{folder}/f{fxx:02d}.png", dpi=90, bbox_inches='tight', facecolor='#020617')
    plt.close()

# --- 1. FIND LATEST *COMPLETED* RUN ---
H_init = None
# Look back 12 hours. We only accept a run if we can find the inventory for hour 18.
for offset in range(12):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        # We test fxx=18. If this exists, the full 1-18 run is likely ready.
        test_H = Herbie(search_time.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=18, save_dir=SAVE_DIR)
        
        if test_H.inventory() is not None:
            H_init = test_H
            print(f"Found latest completed run: {search_time.strftime('%Y-%m-%d %H:00')}Z")
            break
    except:
        continue

if not H_init:
    print("Error: Could not find a completed HRRR run (checked last 12 hours).")
    exit(1)

# --- 2. GENERATION LOOP ---
for fxx in range(1, 19):
    try:
        H = Herbie(H_init.date, model='hrrr', fxx=fxx, save_dir=SAVE_DIR)

        # A. TEMPERATURE
        try:
            ds = safe_get_ds(H, "TMP:2 m")
            data = (ds.t2m - 273.15) * 9/5 + 32
            fig, ax = setup_plot(H, fxx, "Temperature (°F)")
            im = ax.pcolormesh(ds.longitude, ds.latitude, data, transform=ccrs.PlateCarree(), 
                             cmap=ListedColormap(TEMP_COLORS), norm=BoundaryNorm(TEMP_LEVELS, len(TEMP_COLORS)))
            draw_county_labels(ax, ds, data, fmt="{:.0f}")
            # Colorbar with ticks
            cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.9, ticks=TEMP_LEVELS[::2])
            cbar.ax.tick_params(colors='white', labelsize=8)
            save_plot(fig, "frames_temp", fxx)
        except: pass

        # B. WIND CHILL
        try:
            ds = safe_get_ds(H, ":WCHILL:surface")
            var = ds[list(ds.data_vars)[0]]
            data = (var - 273.15) * 9/5 + 32
            fig, ax = setup_plot(H, fxx, "Wind Chill (°F)")
            im = ax.pcolormesh(ds.longitude, ds.latitude, data, transform=ccrs.PlateCarree(), 
                             cmap=ListedColormap(TEMP_COLORS), norm=BoundaryNorm(TEMP_LEVELS, len(TEMP_COLORS)))
            cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.9, ticks=TEMP_LEVELS[::2])
            cbar.ax.tick_params(colors='white', labelsize=8)
            save_plot(fig, "frames_chill", fxx)
        except: pass

        # C. WIND SPEED
        try:
            ds = safe_get_ds(H, ":(UGRD|VGRD):10 m")
            ws = np.sqrt(ds.u10**2 + ds.v10**2) * 2.237
            fig, ax = setup_plot(H, fxx, "Wind Speed (mph)")
            im = ax.pcolormesh(ds.longitude, ds.latitude, ws, transform=ccrs.PlateCarree(), 
                             cmap=ListedColormap(WIND_COLORS), norm=BoundaryNorm(WIND_LEVELS, len(WIND_COLORS)))
            draw_county_labels(ax, ds, ws, fmt="{:.0f}", check_val=10)
            cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.9, ticks=WIND_LEVELS[::2])
            cbar.ax.tick_params(colors='white', labelsize=8)
            save_plot(fig, "frames_wind", fxx)
        except: pass

        # D. WIND GUSTS
        try:
            ds = safe_get_ds(H, ":GUST:surface")
            gusts = ds.gust * 2.237
            fig, ax = setup_plot(H, fxx, "Wind Gusts (mph)")
            im = ax.pcolormesh(ds.longitude, ds.latitude, gusts, transform=ccrs.PlateCarree(), 
                             cmap=ListedColormap(WIND_COLORS), norm=BoundaryNorm(WIND_LEVELS, len(WIND_COLORS)))
            draw_county_labels(ax, ds, gusts, fmt="{:.0f}", check_val=20)
            cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.9, ticks=WIND_LEVELS[::2])
            cbar.ax.tick_params(colors='white', labelsize=8)
            save_plot(fig, "frames_gust", fxx)
        except: pass

        # E. TOTAL PRECIP
        try:
            ds = safe_get_ds(H, ":APCP:surface")
            precip = ds.tp * 0.03937
            fig, ax = setup_plot(H, fxx, "Total Precip (in)")
            im = ax.pcolormesh(ds.longitude, ds.latitude, precip.where(precip >= 0.01), transform=ccrs.PlateCarree(), 
                             cmap=ListedColormap(PRECIP_COLORS), norm=BoundaryNorm(PRECIP_LEVELS, len(PRECIP_COLORS)))
            draw_county_labels(ax, ds, precip, fmt="{:.2f}", check_val=0.05)
            cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.9, ticks=PRECIP_LEVELS)
            cbar.ax.tick_params(colors='white', labelsize=7)
            save_plot(fig, "frames_total_precip", fxx)
        except: pass

        # F. PRECIP TYPE (FIXED MERGE LOGIC)
        try:
            ds = safe_get_ds(H, ":C(RAIN|FREEZ|ICEP|SNOW):surface")
            
            # Ensure variables exist, default to 0 if missing
            crain = ds.crain if 'crain' in ds else xr.zeros_like(ds.longitude)
            cfrzr = ds.cfrzr if 'cfrzr' in ds else xr.zeros_like(ds.longitude)
            icep = ds.icep if 'icep' in ds else xr.zeros_like(ds.longitude)
            # Sometimes snow is 'csnow'
            csnow = ds.csnow if 'csnow' in ds else xr.zeros_like(ds.longitude)

            ptype = crain*1 + cfrzr*2 + icep*3 + csnow*4
            
            fig, ax = setup_plot(H, fxx, "Precipitation Type")
            # Rain (Green), Mix (Orange), Ice (Red), Snow (Midnight Blue)
            cmap = ListedColormap(['#33ff33', '#ff9900', '#ff0000', '#001f3f'])
            norm = BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5], 4)
            im = ax.pcolormesh(ds.longitude, ds.latitude, ptype.where(ptype > 0), transform=ccrs.PlateCarree(), 
                             cmap=cmap, norm=norm)
            
            cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.9, ticks=[1, 2, 3, 4])
            cbar.ax.set_xticklabels(['Rain', 'Mix', 'Ice', 'Snow'], color='white', fontweight='bold')
            save_plot(fig, "frames_precip", fxx)
        except Exception as e: 
            print(f"PType Error f{fxx}: {e}")

        # G. TOTAL SNOWFALL (FIXED MERGE LOGIC)
        try:
            ds = safe_get_ds(H, ":ASNOW:surface")
            snow = ds.asnow * 39.37
            fig, ax = setup_plot(H, fxx, "Total Snowfall (in)")
            im = ax.pcolormesh(ds.longitude, ds.latitude, snow.where(snow >= 0.1), transform=ccrs.PlateCarree(), 
                             cmap=ListedColormap(SNOW_COLORS), norm=BoundaryNorm(SNOW_LEVELS, len(SNOW_COLORS)))
            draw_county_labels(ax, ds, snow, fmt="{:.1f}", check_val=0.5)
            cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.9, ticks=SNOW_LEVELS)
            cbar.ax.tick_params(colors='white', labelsize=7)
            save_plot(fig, "frames_snow", fxx)
        except Exception as e:
            print(f"Snow Error f{fxx}: {e}")

        print(f"Finished Frame f{fxx}")

    except Exception as e:
        print(f"Error on Frame f{fxx}: {e}")
