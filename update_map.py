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

# --- GLOBAL STYLING ---
warnings.filterwarnings("ignore")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Inter', 'Arial', 'Helvetica', 'DejaVu Sans']

# --- CONFIGURATION ---
EXTENTS = [-80.5, -71.5, 38.5, 43.5] 
MAP_CRS = ccrs.LambertConformal(central_longitude=-76.0, central_latitude=41.0)
SAVE_DIR = os.path.join(os.getcwd(), "herbie_data")
SMOOTH_SIGMA = 1.2 

# --- FOLDER SETUP (Run First) ---
folders = ["frames_temp", "frames_chill", "frames_wind", "frames_gust", "frames_snow", "frames_total_precip", "frames_precip"]
for f in folders:
    os.makedirs(f, exist_ok=True)

# --- COLOR PALETTES ---
SNOW_COLORS = ['#081d58', '#253494', '#225ea8', '#1d91c0', '#41b6c4', '#7fcdbb', '#c7e9b4', '#edf8b1', '#fee391', '#fec44f', '#fe9929', '#ec7014', '#cc4c02']
SNOW_LEVELS = [0.1, 1, 2, 3, 4, 6, 8, 12, 18, 24, 30, 36, 48]

WIND_COLORS = ['#1e466e', '#2c69b0', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef', '#e5f5e0', '#a1d99b', '#74c476', '#31a354', '#fd8d3c', '#f03b20', '#bd0026', '#800026']
WIND_LEVELS = [0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100]

PTYPE_COLORS = ['#ffffff00', '#a1d99b', '#41ab5d', '#00441b', '#fbb4b9', '#f768a1', '#7a0177', '#fee5d9', '#ef3b2c', '#67000d', '#d1e5f0', '#4393c3', '#053061']
PTYPE_LEVELS = np.arange(0, 14)

PRECIP_COLORS = ['#ffffff00', '#ffffff', '#99ff33', '#00cc00', '#006600', '#004d66', '#3399ff', '#00ffff', '#9999ff', '#9933ff', '#cc33ff', '#990000', '#cc0000', '#ff3300', '#ff9900', '#cc6600', '#cc9933', '#ffff00', '#ff9999']
PRECIP_LEVELS = [0, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0, 15.0, 20.0, 25.0]

# --- HELPERS ---

def robust_get_data(H, search_list):
    for search in search_list:
        try:
            ds = H.xarray(search, verbose=False)
            if isinstance(ds, list): ds = xr.merge(ds, compat='override')
            if ds is not None: return ds
        except: continue
    return None

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
                    ha='center', va='center', fontsize=6, fontweight='800', color='black', zorder=10)
    except: pass

def setup_plot(H, fxx, datatype):
    fig = plt.figure(figsize=(12, 11), facecolor='#020617')
    ax = plt.axes(projection=MAP_CRS)
    ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='black', linewidth=1.5, zorder=3)
    
    init_dt = H.date.replace(tzinfo=pytz.UTC)
    valid_dt = H.valid_date.replace(tzinfo=pytz.UTC)
    valid_et = valid_dt.astimezone(pytz.timezone('US/Eastern'))
    
    # Clean Header Design
    plt.text(0.5, 1.08, "EPHRATA WEATHER", transform=ax.transAxes, fontsize=22, color='white', fontweight='900', ha='center')
    plt.text(0, 1.02, f"HRRR | {datatype}", transform=ax.transAxes, fontsize=14, color='#94a3b8', ha='left', fontweight='600')
    
    time_str = f"Run: {init_dt.strftime('%m/%d/%Y %H')}Z | F{fxx:02d}\nValid: {valid_dt.strftime('%m/%d/%Y %H')}Z ({valid_et.strftime('%m/%d %I:%M %p')} ET)"
    plt.text(1, 1.02, time_str, transform=ax.transAxes, fontsize=10, color='white', ha='right', weight='bold')
    return fig, ax

# --- FIND LATEST RUN ---
H_init = None
for offset in range(12):
    try:
        t = datetime.utcnow() - timedelta(hours=offset)
        test = Herbie(t.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=18, save_dir=SAVE_DIR)
        if test.inventory() is not None:
            H_init = test; break
    except: continue

if not H_init: exit("No recent HRRR data found.")

# --- MAIN LOOP ---
for fxx in range(1, 19):
    try:
        H = Herbie(H_init.date, model='hrrr', fxx=fxx, save_dir=SAVE_DIR)
        
        # Pull core data for maps and manual calculations
        ds_t = robust_get_data(H, ["TMP:2 m"])
        ds_w = robust_get_data(H, [":(UGRD|VGRD):10 m"])
        
        # 1. TEMPERATURE
        if ds_t:
            temp_f = (ds_t.t2m - 273.15) * 9/5 + 32
            data = gaussian_filter(temp_f, sigma=SMOOTH_SIGMA)
            fig, ax = setup_plot(H, fxx, "Surface Temperature (°F)")
            im = ax.pcolormesh(ds_t.longitude, ds_t.latitude, data, transform=ccrs.PlateCarree(), cmap='magma')
            plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='white')
            draw_labels(ax, ds_t, xr.DataArray(data, coords=ds_t.coords))
            plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='#020617'); plt.close()

        # 2. WIND CHILL (Calculated Fallback)
        if ds_t and ds_w:
            t_f = (ds_t.t2m - 273.15) * 9/5 + 32
            v_mph = np.sqrt(ds_w.u10**2 + ds_w.v10**2) * 2.237
            # NWS Wind Chill Formula
            chill = 35.74 + (0.6215 * t_f) - (35.75 * (v_mph**0.16)) + (0.4275 * t_f * (v_mph**0.16))
            # Wind chill only defined for T <= 50F and V >= 3mph
            chill = xr.where((t_f <= 50) & (v_mph >= 3), chill, t_f)
            data = gaussian_filter(chill, sigma=SMOOTH_SIGMA)
            
            fig, ax = setup_plot(H, fxx, "Wind Chill (°F)")
            im = ax.pcolormesh(ds_t.longitude, ds_t.latitude, data, transform=ccrs.PlateCarree(), cmap='coolwarm')
            plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='white')
            plt.savefig(f"frames_chill/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='#020617'); plt.close()

        # 3. WIND SPEED & GUSTS
        for tag, folder, title in [(":GUST:surface", "frames_gust", "Wind Gusts (mph)"), (":(UGRD|VGRD):10 m", "frames_wind", "Wind Speed (mph)")]:
            ds = robust_get_data(H, [tag])
            if ds:
                val = np.sqrt(ds.u10**2 + ds.v10**2)*2.237 if 'u10' in ds else ds[list(ds.data_vars)[0]]*2.237
                data = gaussian_filter(val, sigma=SMOOTH_SIGMA)
                fig, ax = setup_plot(H, fxx, title)
                im = ax.pcolormesh(ds.longitude, ds.latitude, data, transform=ccrs.PlateCarree(), cmap=ListedColormap(WIND_COLORS), norm=BoundaryNorm(WIND_LEVELS, 14))
                plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='white')
                draw_labels(ax, ds, xr.DataArray(data, coords=ds.coords), check_val=15)
                plt.savefig(f"{folder}/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='#020617'); plt.close()

        # 4. TOTAL PRECIP 
        ds_tp = robust_get_data(H, [":APCP:surface"])
        if ds_tp:
            data = gaussian_filter(ds_tp.tp * 0.03937, sigma=SMOOTH_SIGMA)
            fig, ax = setup_plot(H, fxx, "Total Precipitation (in)")
            im = ax.pcolormesh(ds_tp.longitude, ds_tp.latitude, np.where(data > 0.01, data, np.nan), 
                               transform=ccrs.PlateCarree(), cmap=ListedColormap(PRECIP_COLORS), norm=BoundaryNorm(PRECIP_LEVELS, len(PRECIP_COLORS)))
            plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='white')
            draw_labels(ax, ds_tp, xr.DataArray(data, coords=ds_tp.coords), fmt="{:.2f}", check_val=0.1)
            plt.savefig(f"frames_total_precip/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='#020617'); plt.close()

        # 5. PRECIP TYPE (Intensity Based)
        ds_pt = robust_get_data(H, [":(REFC|CRAIN|CSNOW|CFRZR|CICEP):"])
        if ds_pt:
            refc = gaussian_filter(ds_pt.refc.values, sigma=SMOOTH_SIGMA)
            intensity = np.clip((refc / 15).astype(int), 0, 2)
            final_map = np.zeros_like(refc)
            for i, vname in enumerate(['crain', 'cfrzr', 'icep', 'csnow'], 0):
                if vname in ds_pt: final_map = np.where(ds_pt[vname] == 1, (i*3)+1 + intensity, final_map)
            fig, ax = setup_plot(H, fxx, "Precip Type & Intensity")
            im = ax.pcolormesh(ds_pt.longitude, ds_pt.latitude, np.ma.masked_where(final_map==0, final_map),
                               transform=ccrs.PlateCarree(), cmap=ListedColormap(PTYPE_COLORS), norm=BoundaryNorm(PTYPE_LEVELS, 14))
            cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8, ticks=[2, 5, 8, 11])
            cbar.ax.set_xticklabels(['Rain', 'Mix', 'Ice', 'Snow'], color='white', fontweight='bold')
            plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='#020617'); plt.close()

        # 6. SNOWFALL
        ds_s = robust_get_data(H, [":ASNOW:surface", ":WEASD:surface", ":SNOD:surface"])
        if ds_s:
            data = gaussian_filter(ds_s[list(ds_s.data_vars)[0]] * 39.37, sigma=SMOOTH_SIGMA)
            fig, ax = setup_plot(H, fxx, "Total Snowfall (in)")
            im = ax.pcolormesh(ds_s.longitude, ds_s.latitude, np.where(data > 0.1, data, np.nan), 
                               transform=ccrs.PlateCarree(), cmap=ListedColormap(SNOW_COLORS), norm=BoundaryNorm(SNOW_LEVELS, 13))
            plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='white')
            draw_labels(ax, ds_s, xr.DataArray(data, coords=ds_s.coords), fmt="{:.1f}", check_val=0.1)
            plt.savefig(f"frames_snow/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='#020617'); plt.close()

        print(f"Finished Frame f{fxx}")
    except Exception as e:
        print(f"Error on f{fxx}: {e}")
