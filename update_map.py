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
SMOOTH_SIGMA = 1.2  # Smoothing intensity for visual appeal

# --- COLOR PALETTES ---
# Intensity-based P-Type (Tropical Tidbits Inspired)
PTYPE_COLORS = [
    '#ffffff00', # Transparent
    '#a1d99b', '#41ab5d', '#00441b', # Rain (L, M, H)
    '#fbb4b9', '#f768a1', '#7a0177', # Mix/Frz (L, M, H)
    '#fee5d9', '#ef3b2c', '#67000d', # Ice/Sleet (L, M, H)
    '#d1e5f0', '#4393c3', '#053061'  # Snow (L, M, H)
]
PTYPE_LEVELS = np.arange(0, 14)

# High-Visibility Snow Scale (Dark Blues to Reds)
SNOW_COLORS = ['#08306b', '#08519c', '#2171b5', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef', '#deebf7', '#ffffcc', '#ffeda0', '#feb24c', '#f03b20', '#bd0026']
SNOW_LEVELS = [0.1, 1, 2, 3, 4, 6, 8, 12, 18, 24, 30, 36, 48]

WIND_COLORS = ['#1e466e', '#2c69b0', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef', '#e5f5e0', '#a1d99b', '#74c476', '#31a354', '#fd8d3c', '#f03b20', '#bd0026', '#800026']
WIND_LEVELS = [0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100]

# --- HELPERS ---

def safe_get_ds(H, search_str):
    ds = H.xarray(search_str, verbose=False)
    if isinstance(ds, list):
        return xr.merge(ds, compat='override')
    return ds

def draw_labels(ax, ds, data_var, fmt="{:.0f}", check_val=None):
    try:
        reader = shpreader.Reader(shpreader.natural_earth(resolution='10m', category='cultural', name='admin_2_counties'))
        counties = list(reader.records())
        for county in counties:
            lon, lat = county.geometry.centroid.x, county.geometry.centroid.y
            if lon < EXTENTS[0] or lon > EXTENTS[1] or lat < EXTENTS[2] or lat > EXTENTS[3]: continue
            val = float(data_var.sel(latitude=lat, longitude=lon, method='nearest').values)
            if check_val is not None and val < check_val: continue
            ax.text(lon, lat, fmt.format(val), transform=ccrs.PlateCarree(),
                    ha='center', va='center', fontsize=6, fontweight='bold', color='black', zorder=10)
    except: pass

def setup_plot(H, fxx, datatype):
    fig = plt.figure(figsize=(12, 11), facecolor='#020617')
    ax = plt.axes(projection=MAP_CRS)
    ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='black', linewidth=1.5, zorder=3)
    
    init_dt = H.date.replace(tzinfo=pytz.UTC)
    valid_dt = H.valid_date.replace(tzinfo=pytz.UTC)
    valid_et = valid_dt.astimezone(pytz.timezone('US/Eastern'))
    
    # Header Branding
    plt.text(0.5, 1.06, "EPHRATA WEATHER", transform=ax.transAxes, fontsize=18, color='white', fontweight='black', ha='center')
    plt.text(0, 1.02, f"HRRR: {datatype}", transform=ax.transAxes, fontsize=14, color='white', ha='left', fontweight='bold')
    
    time_str = f"Run: {init_dt.strftime('%m/%d/%Y %H')}Z | F{fxx:02d}\nValid: {valid_dt.strftime('%m/%d/%Y %H')}Z ({valid_et.strftime('%m/%d %I:%M %p')} ET)"
    plt.text(1, 1.02, time_str, transform=ax.transAxes, fontsize=9, color='white', ha='right', weight='bold')
    return fig, ax

# --- INITIALIZATION ---
folders = ["frames_temp", "frames_chill", "frames_wind", "frames_gust", "frames_snow", "frames_total_precip", "frames_precip"]
for f in folders: os.makedirs(f, exist_ok=True)

H_init = None
for offset in range(12):
    try:
        t = datetime.utcnow() - timedelta(hours=offset)
        test = Herbie(t.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=18, save_dir=SAVE_DIR)
        if test.inventory() is not None:
            H_init = test; break
    except: continue

# --- MAIN LOOP ---
for fxx in range(1, 19):
    try:
        H = Herbie(H_init.date, model='hrrr', fxx=fxx, save_dir=SAVE_DIR)
        
        # 1. TEMPERATURE
        try:
            ds = safe_get_ds(H, "TMP:2 m")
            temp = gaussian_filter((ds.t2m - 273.15) * 9/5 + 32, sigma=SMOOTH_SIGMA)
            fig, ax = setup_plot(H, fxx, "2m Temperature (°F)")
            im = ax.pcolormesh(ds.longitude, ds.latitude, temp, transform=ccrs.PlateCarree(), cmap='magma')
            plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='white')
            draw_labels(ax, ds, xr.DataArray(temp, coords=ds.coords))
            plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='#020617'); plt.close()
        except Exception as e: print(f"Temp Error: {e}")

        # 2. PRECIP TYPE (Intensity Based)
        try:
            ds = safe_get_ds(H, ":(REFC|CRAIN|CSNOW|CFRZR|CICEP):")
            refc = gaussian_filter(ds.refc.values, sigma=SMOOTH_SIGMA)
            intensity = np.clip((refc / 15).astype(int), 0, 2)
            final_map = np.zeros_like(refc)
            if 'crain' in ds: final_map = np.where(ds.crain == 1, 1 + intensity, final_map)
            if 'cfrzr' in ds: final_map = np.where(ds.cfrzr == 1, 4 + intensity, final_map)
            if 'icep' in ds:  final_map = np.where(ds.icep == 1, 7 + intensity, final_map)
            if 'csnow' in ds: final_map = np.where(ds.csnow == 1, 10 + intensity, final_map)
            fig, ax = setup_plot(H, fxx, "Precipitation Type & Intensity")
            im = ax.pcolormesh(ds.longitude, ds.latitude, np.ma.masked_where(final_map==0, final_map),
                               transform=ccrs.PlateCarree(), cmap=ListedColormap(PTYPE_COLORS), norm=BoundaryNorm(PTYPE_LEVELS, 14))
            cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8, ticks=[2, 5, 8, 11])
            cbar.ax.set_xticklabels(['Rain', 'Mix', 'Ice', 'Snow'], color='white', fontweight='bold')
            plt.text(0.99, -0.04, "Colors: TropicalTidbits.com Style", transform=ax.transAxes, color='gray', fontsize=7, ha='right')
            plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='#020617'); plt.close()
        except: pass

        # 3. WIND SPEED & GUSTS
        for tag, folder, title in [(":GUST:", "frames_gust", "Wind Gusts (mph)"), (":(UGRD|VGRD):10 m", "frames_wind", "Wind Speed (mph)")]:
            try:
                ds = safe_get_ds(H, tag)
                val = np.sqrt(ds.u10**2 + ds.v10**2)*2.237 if 'u10' in ds else ds[list(ds.data_vars)[0]]*2.237
                data = gaussian_filter(val, sigma=SMOOTH_SIGMA)
                fig, ax = setup_plot(H, fxx, title)
                im = ax.pcolormesh(ds.longitude, ds.latitude, data, transform=ccrs.PlateCarree(), cmap=ListedColormap(WIND_COLORS), norm=BoundaryNorm(WIND_LEVELS, 14))
                plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='white')
                draw_labels(ax, ds, xr.DataArray(data, coords=ds.coords), check_val=15)
                plt.savefig(f"{folder}/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='#020617'); plt.close()
            except: continue

        # 4. TOTAL SNOWFALL
        try:
            ds = safe_get_ds(H, ":(ASNOW|WEASD|SNOD):surface")
            snow = gaussian_filter(ds[list(ds.data_vars)[0]] * 39.37, sigma=SMOOTH_SIGMA)
            fig, ax = setup_plot(H, fxx, "Total Snowfall (in)")
            im = ax.pcolormesh(ds.longitude, ds.latitude, np.where(snow > 0.1, snow, np.nan), 
                               transform=ccrs.PlateCarree(), cmap=ListedColormap(SNOW_COLORS), norm=BoundaryNorm(SNOW_LEVELS, 13))
            plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='white')
            draw_labels(ax, ds, xr.DataArray(snow, coords=ds.coords), fmt="{:.1f}", check_val=0.1)
            plt.savefig(f"frames_snow/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='#020617'); plt.close()
        except: pass

        print(f"Frame f{fxx} Ready")
    except: continue
