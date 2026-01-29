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
from scipy.ndimage import gaussian_filter
import warnings

warnings.filterwarnings("ignore")

# --- CONFIGURATION ---
EXTENTS = [-80.5, -71.5, 38.5, 43.5] 
MAP_CRS = ccrs.LambertConformal(central_longitude=-76.0, central_latitude=41.0)
SAVE_DIR = os.path.join(os.getcwd(), "herbie_data")
SMOOTH_SIGMA = 1.0  # Increase for more smoothing, 0 for raw data

# --- COLOR PALETTES ---
# Tropical Tidbits inspired P-Type (Rain=Green, Mix=Pink, Ice=Red, Snow=Blue)
TT_PTYPE_COLORS = ['#00ad00', '#ff00f3', '#ff0000', '#0000ff'] 
TT_PTYPE_LEVELS = [0.5, 1.5, 2.5, 3.5, 4.5]

WIND_COLORS = ['#1e466e', '#2c69b0', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef', '#e5f5e0', '#a1d99b', '#74c476', '#31a354', '#fd8d3c', '#f03b20', '#bd0026', '#800026', '#49006a']
WIND_LEVELS = [0, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 140, 160]

# --- HELPERS ---

def safe_get_ds(H, search_str):
    ds = H.xarray(search_str, verbose=False)
    if isinstance(ds, list):
        return xr.merge(ds, compat='override')
    return ds

def smooth_data(data):
    """Applies Gaussian smoothing to make the data more visually appealing."""
    return gaussian_filter(data, sigma=SMOOTH_SIGMA)

def setup_plot(H, fxx, datatype):
    fig = plt.figure(figsize=(12, 11), facecolor='#020617')
    ax = plt.axes(projection=MAP_CRS)
    ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
    
    # States and Borders
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='black', linewidth=1.5, zorder=3)
    
    init_dt = H.date.replace(tzinfo=pytz.UTC)
    valid_dt = H.valid_date.replace(tzinfo=pytz.UTC)
    valid_et = valid_dt.astimezone(pytz.timezone('US/Eastern'))
    
    # Brand Bar at top to avoid legend overlap
    plt.text(0.5, 1.06, "EPHRATA WEATHER", transform=ax.transAxes, 
             fontsize=18, color='white', fontweight='black', ha='center')
    
    # Title Layout
    plt.text(0, 1.02, f"HRRR: {datatype}", transform=ax.transAxes, 
             fontsize=14, color='white', fontweight='bold', ha='left')
    
    run_info = f"Run: {init_dt.strftime('%m/%d/%Y %H')}Z | F{fxx:02d}\n"
    valid_info = f"Valid: {valid_dt.strftime('%m/%d/%Y %H')}Z ({valid_et.strftime('%m/%d %I:%M %p')} ET)"
    plt.text(1, 1.02, run_info + valid_info, transform=ax.transAxes, 
             fontsize=9, color='white', ha='right', weight='bold')

    return fig, ax

def save_plot(fig, folder, fxx):
    os.makedirs(folder, exist_ok=True)
    plt.savefig(f"{folder}/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='#020617')
    plt.close()

# --- 1. SETUP FOLDERS ---
folders = ["frames_temp", "frames_chill", "frames_wind", "frames_gust", "frames_snow", "frames_total_precip", "frames_precip"]
for f in folders: os.makedirs(f, exist_ok=True)

# --- 2. FIND LATEST COMPLETE RUN ---
H_init = None
for offset in range(12):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        test_H = Herbie(search_time.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=18, save_dir=SAVE_DIR)
        if test_H.inventory() is not None:
            H_init = test_H
            print(f"Loading Run: {search_time.strftime('%Y-%m-%d %H:00')}Z")
            break
    except: continue

if not H_init: exit("No complete data found.")

# --- 3. GENERATION LOOP ---
for fxx in range(1, 19):
    try:
        H = Herbie(H_init.date, model='hrrr', fxx=fxx, save_dir=SAVE_DIR)

        # A. PRECIP TYPE (Tropical Tidbits Style)
        try:
            ds = safe_get_ds(H, ":C(RAIN|FREEZ|ICEP|SNOW):surface")
            crain = ds.crain if 'crain' in ds else xr.zeros_like(ds.longitude)
            cfrzr = ds.cfrzr if 'cfrzr' in ds else xr.zeros_like(ds.longitude)
            icep = ds.icep if 'icep' in ds else xr.zeros_like(ds.longitude)
            csnow = ds.csnow if 'csnow' in ds else xr.zeros_like(ds.longitude)
            
            ptype = crain*1 + cfrzr*2 + icep*3 + csnow*4
            fig, ax = setup_plot(H, fxx, "Precipitation Type")
            
            im = ax.pcolormesh(ds.longitude, ds.latitude, ptype.where(ptype > 0), transform=ccrs.PlateCarree(), 
                             cmap=ListedColormap(TT_PTYPE_COLORS), norm=BoundaryNorm(TT_PTYPE_LEVELS, 4))
            
            cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8, ticks=[1, 2, 3, 4])
            cbar.ax.set_xticklabels(['Rain', 'Mix/Frz Rain', 'Sleet/Ice', 'Snow'], color='white', fontweight='bold')
            plt.text(0.99, -0.04, "Colors: TropicalTidbits.com", transform=ax.transAxes, color='gray', fontsize=7, ha='right')
            save_plot(fig, "frames_precip", fxx)
        except Exception as e: print(f"PType Error: {e}")

        # B. SNOWFALL (Fixing 'asnow' attribute error)
        try:
            # Fallback logic for snowfall variables
            snow_vars = [":ASNOW:surface", ":WEASD:surface", ":SNOD:surface"]
            ds = None
            for sv in snow_vars:
                try:
                    ds = safe_get_ds(H, sv)
                    break
                except: continue
            
            if ds:
                var_name = list(ds.data_vars)[0]
                snow_data = ds[var_name] * 39.37 # Meters to Inches
                # Apply smoothing
                display_snow = smooth_data(snow_data)
                
                fig, ax = setup_plot(H, fxx, "Total Snowfall (in)")
                im = ax.pcolormesh(ds.longitude, ds.latitude, np.where(display_snow > 0.1, display_snow, np.nan), 
                                 transform=ccrs.PlateCarree(), cmap='Blues', vmin=0, vmax=12)
                plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='white')
                save_plot(fig, "frames_snow", fxx)
        except Exception as e: print(f"Snow Error: {e}")

        # C. WIND CHILL
        try:
            ds = safe_get_ds(H, ":WCHILL:surface")
            var = ds[list(ds.data_vars)[0]]
            data = smooth_data((var - 273.15) * 9/5 + 32)
            fig, ax = setup_plot(H, fxx, "Wind Chill (°F)")
            im = ax.pcolormesh(ds.longitude, ds.latitude, data, transform=ccrs.PlateCarree(), cmap='coolwarm')
            plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='white')
            save_plot(fig, "frames_chill", fxx)
        except: pass

        print(f"Finished Frame f{fxx}")

    except Exception as e:
        print(f"Fatal error f{fxx}: {e}")
