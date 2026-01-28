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
EXTENTS = [-80.5, -71.5, 38.5, 43.5] # Focused Northeast (PA/NY/NJ/MD/DE)
MAP_CRS = ccrs.LambertConformal(central_longitude=-76.0, central_latitude=41.0)

# Pre-load County Shapes for labels and borders
try:
    reader = shpreader.Reader(shpreader.natural_earth(resolution='10m', category='cultural', name='admin_2_counties'))
    counties = list(reader.records())
except:
    counties = []

# --- HELPERS ---

def get_data_array(ds):
    """Safely extracts the first data variable regardless of its name to prevent crashes."""
    if isinstance(ds, list):
        ds = xr.merge(ds, compat='override')
    # Filter out coordinate variables to find the actual data
    vars = [v for v in ds.data_vars]
    if not vars:
        raise ValueError("No data variables found in dataset")
    return ds, ds[vars[0]]

def draw_county_labels(ax, ds, data_var, fmt="{:.0f}", color='black', check_val=None):
    """Samples data at the center of each county and draws a label."""
    view_w, view_e, view_s, view_n = EXTENTS
    for county in counties:
        geom = county.geometry
        bounds = geom.bounds
        # Only process counties in the current view
        if bounds[2] < view_w or bounds[0] > view_e or bounds[3] < view_s or bounds[1] > view_n:
            continue
        
        centroid = geom.centroid
        lon, lat = centroid.x, centroid.y
        
        try:
            # Sample the data at the county center
            val = float(data_var.sel(latitude=lat, longitude=lon, method='nearest').values)
            if check_val is not None and val < check_val: 
                continue
            
            ax.text(lon, lat, fmt.format(val), transform=ccrs.PlateCarree(),
                    ha='center', va='center', fontsize=6, fontweight='bold',
                    color=color, clip_on=True, zorder=10)
        except:
            continue

def add_precip_legend(ax):
    """Adds a legend identifying Rain, Snow, Sleet, and Freezing Rain."""
    handles = [
        mpatches.Patch(color='#e67e22', label='Rain'),
        mpatches.Patch(color='#08306b', label='Snow'),
        mpatches.Patch(color='#6a51a3', label='Sleet'),
        mpatches.Patch(color='#ae017e', label='Frz. Rain')
    ]
    ax.legend(handles=handles, loc='lower left', fontsize=8, facecolor='white', framealpha=0.9)

def get_precip_cmap():
    """Flipped snow intensity: Light blue = Light Snow, Dark blue = Heavy Snow."""
    colors = ['#ffffff00'] # Transparent for no precip
    colors.extend(['#3498db', '#2ecc71', '#f1c40f', '#e67e22', '#e74c3c']) # Rain
    colors.extend(['#deebf7', '#6baed6', '#2171b5', '#08519c', '#08306b']) # SNOW (Light to Dark)
    colors.extend(['#efedf5', '#9e9ac8', '#6a51a3', '#54278f', '#3f007d']) # Sleet
    colors.extend(['#fbb4b9', '#f768a1', '#dd3497', '#ae017e', '#7a0177']) # Ice
    return ListedColormap(colors)

def setup_map(title):
    fig = plt.figure(figsize=(12, 10))
    ax = plt.axes(projection=MAP_CRS)
    ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
    
    # Add States and County borders
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='black', linewidth=1.2, zorder=2)
    ax.add_feature(cfeature.ShapelyFeature([c.geometry for c in counties], ccrs.PlateCarree()), 
                   facecolor='none', edgecolor='gray', linewidth=0.4, alpha=0.6, zorder=1)
    
    plt.title(title, loc='left', fontweight='bold', fontsize=12)
    return fig, ax

# --- DATA SEARCH ---
H_init = None
print("Searching for latest HRRR run...")
for offset in range(24):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        test_H = Herbie(search_time.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=1)
        if test_H.inventory() is not None:
            H_init = test_H
            print(f"Using run: {search_time.strftime('%Y-%m-%d %H:00')}")
            break
    except: continue

if not H_init: 
    print("Could not find any HRRR data!")
    exit(1)

# Ensure folders exist
folders = ["frames_temp", "frames_precip", "frames_wind", "frames_gust", "frames_snow", "frames_total_precip"]
for folder in folders: os.makedirs(folder, exist_ok=True)

# --- GENERATION LOOP ---
for fxx in range(1, 19):
    try:
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=fxx)
        utc_v = H.valid_date.replace(tzinfo=pytz.UTC)
        et_v = utc_v.astimezone(pytz.timezone('US/Eastern'))
        t_str = f"Valid: {et_v.strftime('%m/%d %I:%M %p ET')}"

        # 1. TEMPERATURE (2m)
        try:
            ds_raw = H.xarray("TMP:2 m", verbose=False)
            ds, temp_k = get_data_array(ds_raw)
            temp_f = (temp_k - 273.15) * 9/5 + 32
            fig, ax = setup_map(f"HRRR 2m Temperature (°F) | {t_str}")
            im = ax.pcolormesh(ds.longitude, ds.latitude, temp_f, transform=ccrs.PlateCarree(), cmap='jet', vmin=0, vmax=100)
            plt.colorbar(im, fraction=0.046, pad=0.04)
            draw_county_labels(ax, ds, temp_f)
            plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Temp Error: {e}")

        # 2. PRECIP TYPE (Flipped Snow Intensity)
        try:
            ds_p = H.xarray(":(REFC|CRAIN|CSNOW|CFRZR|CICEP):", verbose=False)
            if isinstance(ds_p, list): ds_p = xr.merge(ds_p)
            refc = ds_p.refc.values
            intensity = np.clip((refc / 10).astype(int), 0, 4)
            final_map = np.zeros_like(refc)
            if 'crain' in ds_p: final_map = np.where(ds_p.crain == 1, 1 + intensity, final_map)
            if 'csnow' in ds_p: final_map = np.where(ds_p.csnow == 1, 6 + intensity, final_map)
            if 'cicep' in ds_p: final_map = np.where(ds_p.cicep == 1, 11 + intensity, final_map)
            if 'cfrzr' in ds_p: final_map = np.where(ds_p.cfrzr == 1, 16 + intensity, final_map)
            fig, ax = setup_map(f"HRRR Precip Type | {t_str}")
            cmap = get_precip_cmap()
            ax.pcolormesh(ds_p.longitude, ds_p.latitude, np.ma.masked_where(final_map==0, final_map),
                          transform=ccrs.PlateCarree(), cmap=cmap, norm=BoundaryNorm(np.arange(0, 22), cmap.N))
            add_precip_legend(ax)
            plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"P-Type Error: {e}")

        # 3. WIND SPEED (10m)
        try:
            ds_w = H.xarray(":(UGRD|VGRD):10 m", verbose=False)
            if isinstance(ds_w, list): ds_w = xr.merge(ds_w)
            u = ds_w['u10'] if 'u10' in ds_w else ds_w['u']
            v = ds_w['v10'] if 'v10' in ds_w else ds_w['v']
            ws = np.sqrt(u**2 + v**2) * 1.944
            fig, ax = setup_map(f"HRRR Wind Speed (kts) | {t_str}")
            im = ax.pcolormesh(ds_w.longitude, ds_w.latitude, ws, transform=ccrs.PlateCarree(), cmap='BuPu', vmin=0, vmax=50)
            plt.colorbar(im, fraction=0.046, pad=0.04)
            draw_county_labels(ax, ds_w, ws, check_val=5)
            plt.savefig(f"frames_wind/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Wind Error: {e}")

        # 4. WIND GUSTS
        try:
            ds_raw = H.xarray(":GUST:surface", verbose=False)
            ds, gust_var = get_data_array(ds_raw)
            gust_kts = gust_var * 1.944
            fig, ax = setup_map(f"HRRR Wind Gust (kts) | {t_str}")
            im = ax.pcolormesh(ds.longitude, ds.latitude, gust_kts, transform=ccrs.PlateCarree(), cmap='Reds', vmin=0, vmax=60)
            plt.colorbar(im, fraction=0.046, pad=0.04)
            draw_county_labels(ax, ds, gust_kts, check_val=15)
            plt.savefig(f"frames_gust/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Gust Error: {e}")

        # 5. TOTAL SNOWFALL (NWS Scale)
        try:
            ds_raw = H.xarray(":ASNOW:surface", verbose=False)
            ds, snow_var = get_data_array(ds_raw)
            snow_in = snow_var * 39.37
            snow_colors = ['#ffffff00', '#d1e3f3', '#95cbee', '#529dcc', '#1c64a5', '#08306b', '#ffffcc', '#ffeda0', '#feb24c', '#f03b20', '#bd0026', '#7f0000', '#4d0000', '#cbc9e2', '#9e9ac8', '#6a51a3']
            snow_levels = [0, 0.1, 1, 2, 3, 4, 6, 8, 12, 18, 24, 30, 36, 48, 60, 72, 100]
            snow_cmap = ListedColormap(snow_colors)
            fig, ax = setup_map(f"HRRR Total Snow (in) | {t_str}")
            im = ax.pcolormesh(ds.longitude, ds.latitude, snow_in, transform=ccrs.PlateCarree(), cmap=snow_cmap, norm=BoundaryNorm(snow_levels, len(snow_colors)))
            plt.colorbar(im, fraction=0.046, pad=0.04)
            draw_county_labels(ax, ds, snow_in, fmt="{:.1f}", check_val=0.1)
            plt.savefig(f"frames_snow/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Snow Error: {e}")

        # 6. TOTAL PRECIP (High-Res Scale)
        try:
            ds_raw = H.xarray(":APCP:surface", verbose=False)
            ds, precip_var = get_data_array(ds_raw)
            precip_in = precip_var * 0.03937
            tp_colors = ['#ffffff00', '#ffffff', '#99ff33', '#00cc00', '#006600', '#004d66', '#3399ff', '#00ffff', '#9999ff', '#9933ff', '#cc33ff', '#990000', '#cc0000', '#ff3300', '#ff9900', '#cc6600', '#cc9933', '#ffff00', '#ff9999']
            tp_levels = [0, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 50.0]
            fig, ax = setup_map(f"HRRR Total Precip (in) | {t_str}")
            im = ax.pcolormesh(ds.longitude, ds.latitude, precip_in, transform=ccrs.PlateCarree(), cmap=ListedColormap(tp_colors), norm=BoundaryNorm(tp_levels, len(tp_colors)))
            plt.colorbar(im, fraction=0.046, pad=0.04)
            draw_county_labels(ax, ds, precip_in, fmt="{:.2f}", check_val=0.01)
            plt.savefig(f"frames_total_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Total Precip Error: {e}")

        print(f"Completed f{fxx}")
    except Exception as e:
        print(f"Fatal error on f{fxx}: {e}")
