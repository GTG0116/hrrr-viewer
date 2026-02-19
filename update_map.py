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
TEMP_COLORS = [
    '#3b0764', '#5b21b6', '#1e3a8a', '#1d4ed8', '#0ea5e9',
    '#06b6d4', '#10b981', '#22c55e', '#84cc16', '#eab308',
    '#f97316', '#ef4444', '#dc2626', '#991b1b', '#7f1d1d']
TEMP_LEVELS = [-30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]

CHILL_COLORS = [
    '#312e81', '#3730a3', '#4338ca', '#1d4ed8', '#2563eb',
    '#3b82f6', '#06b6d4', '#14b8a6', '#22c55e', '#84cc16',
    '#eab308', '#f97316', '#ef4444', '#dc2626', '#991b1b']
CHILL_LEVELS = [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

SNOW_COLORS = [
    '#bae6fd', '#7dd3fc', '#38bdf8', '#0ea5e9', '#0284c7',
    '#0369a1', '#1e40af', '#a855f7', '#7c3aed', '#f97316',
    '#ef4444', '#dc2626', '#991b1b']
SNOW_LEVELS = [0.1, 1, 2, 3, 4, 6, 8, 12, 18, 24, 30, 36, 48]

WIND_COLORS = [
    '#dbeafe', '#93c5fd', '#60a5fa', '#3b82f6', '#2563eb',
    '#1d4ed8', '#10b981', '#059669', '#eab308', '#f97316',
    '#ef4444', '#dc2626', '#991b1b', '#7f1d1d']
WIND_LEVELS = [0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100]

PTYPE_COLORS = [
    '#ffffff00',                                        # 0:    transparent (no precipitation)
    '#86efac', '#22c55e', '#eab308', '#f97316', '#dc2626',      # 1-5:   Rain  (lt green → green → yellow → orange → red)
    '#fda4af', '#fb923c', '#ef4444', '#991b1b', '#4c0519',      # 6-10:  Frzg Rain (lt pink → salmon → red → dark red → maroon)
    '#e9d5ff', '#c084fc', '#9333ea', '#6b21a8', '#3b0764',      # 11-15: Ice   (lt lavender → lavender → purple → dk purple → deepest)
    '#bae6fd', '#38bdf8', '#0d9488', '#9333ea', '#f9a8d4',      # 16-20: Snow  (lt sky → med blue → teal → purple → pink)
]
PTYPE_LEVELS = np.arange(0, 22)

PRECIP_COLORS = [
    '#ffffff00', '#e0f2fe', '#7dd3fc', '#38bdf8', '#0ea5e9',
    '#0284c7', '#0369a1', '#1d4ed8', '#a855f7', '#7c3aed',
    '#6d28d9', '#ef4444', '#dc2626', '#b91c1c', '#f97316',
    '#eab308', '#84cc16', '#fbbf24', '#f472b6']
PRECIP_LEVELS = [0, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0, 15.0, 20.0, 25.0]

# --- CITY COORDINATES (lon, lat) ---
CITIES = {
    'Philadelphia': (-75.165, 39.953),
    'Pittsburgh':   (-79.996, 40.441),
    'Buffalo':      (-78.878, 42.886),
    'NYC':          (-74.006, 40.713),
    'Lancaster':    (-76.306, 40.038),
    'Harrisburg':   (-76.887, 40.273),
    'Allentown':    (-75.471, 40.602),
    'Scranton':     (-75.662, 41.409),
    'Erie':         (-80.085, 42.129),
    'Baltimore':    (-76.612, 39.290),
}

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

def draw_city_labels(ax, data_array=None, fmt="{:.0f}", suffix=""):
    """Plot city name markers; optionally annotate with the nearest-grid-point value."""
    for city, (lon, lat) in CITIES.items():
        if lon < EXTENTS[0] or lon > EXTENTS[1] or lat < EXTENTS[2] or lat > EXTENTS[3]:
            continue
        ax.plot(lon, lat, 'o', markersize=3.5, color='#0f172a', markeredgewidth=0.5,
                transform=ccrs.PlateCarree(), zorder=12)
        label = city
        if data_array is not None:
            try:
                val = float(data_array.sel(latitude=lat, longitude=lon, method='nearest').values)
                label = f"{city}\n{fmt.format(val)}{suffix}"
            except:
                pass
        ax.text(lon + 0.13, lat + 0.1, label, transform=ccrs.PlateCarree(),
                ha='left', va='bottom', fontsize=6.5, fontweight='bold', color='#0f172a', zorder=13,
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.72, edgecolor='none'))

def setup_plot(H, fxx, datatype):
    fig = plt.figure(figsize=(12, 11), facecolor='white')
    ax = plt.axes(projection=MAP_CRS)
    ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_2_counties_lakes', '10m',
                   facecolor='none', edgecolor='#94a3b8', linewidth=0.3), zorder=2)
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='#334155', linewidth=1.5, zorder=3)

    init_dt = H.date.replace(tzinfo=pytz.UTC)
    valid_dt = H.valid_date.replace(tzinfo=pytz.UTC)
    valid_et = valid_dt.astimezone(pytz.timezone('US/Eastern'))

    # Clean Header Design
    plt.text(0.5, 1.08, "EPHRATA WEATHER", transform=ax.transAxes, fontsize=22, color='#0f172a', fontweight='900', ha='center')
    plt.text(0, 1.02, f"HRRR | {datatype}", transform=ax.transAxes, fontsize=14, color='#475569', ha='left', fontweight='600')

    time_str = f"Run: {init_dt.strftime('%m/%d/%Y %H')}Z | F{fxx:02d}\nValid: {valid_dt.strftime('%m/%d/%Y %H')}Z ({valid_et.strftime('%m/%d %I:%M %p')} ET)"
    plt.text(1, 1.02, time_str, transform=ax.transAxes, fontsize=10, color='#334155', ha='right', weight='bold')
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
            fig, ax = setup_plot(H, fxx, "Surface Temperature (\u00b0F)")
            im = ax.pcolormesh(ds_t.longitude, ds_t.latitude, data, transform=ccrs.PlateCarree(),
                               cmap=ListedColormap(TEMP_COLORS), norm=BoundaryNorm(TEMP_LEVELS, len(TEMP_COLORS)))
            plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='#334155')
            draw_labels(ax, ds_t, xr.DataArray(data, coords=ds_t.coords))
            draw_city_labels(ax, xr.DataArray(data, coords=ds_t.coords), fmt="{:.0f}", suffix="°")
            plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='white'); plt.close()

        # 2. WIND CHILL (Calculated Fallback)
        if ds_t and ds_w:
            t_f = (ds_t.t2m - 273.15) * 9/5 + 32
            v_mph = np.sqrt(ds_w.u10**2 + ds_w.v10**2) * 2.237
            # NWS Wind Chill Formula
            chill = 35.74 + (0.6215 * t_f) - (35.75 * (v_mph**0.16)) + (0.4275 * t_f * (v_mph**0.16))
            # Wind chill only defined for T <= 50F and V >= 3mph
            chill = xr.where((t_f <= 50) & (v_mph >= 3), chill, t_f)
            data = gaussian_filter(chill, sigma=SMOOTH_SIGMA)

            fig, ax = setup_plot(H, fxx, "Wind Chill (\u00b0F)")
            im = ax.pcolormesh(ds_t.longitude, ds_t.latitude, data, transform=ccrs.PlateCarree(),
                               cmap=ListedColormap(CHILL_COLORS), norm=BoundaryNorm(CHILL_LEVELS, len(CHILL_COLORS)))
            plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='#334155')
            plt.savefig(f"frames_chill/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='white'); plt.close()

        # 3. WIND SPEED & GUSTS
        for tag, folder, title in [(":GUST:surface", "frames_gust", "Wind Gusts (mph)"), (":(UGRD|VGRD):10 m", "frames_wind", "Wind Speed (mph)")]:
            ds = robust_get_data(H, [tag])
            if ds:
                val = np.sqrt(ds.u10**2 + ds.v10**2)*2.237 if 'u10' in ds else ds[list(ds.data_vars)[0]]*2.237
                data = gaussian_filter(val, sigma=SMOOTH_SIGMA)
                fig, ax = setup_plot(H, fxx, title)
                im = ax.pcolormesh(ds.longitude, ds.latitude, data, transform=ccrs.PlateCarree(), cmap=ListedColormap(WIND_COLORS), norm=BoundaryNorm(WIND_LEVELS, 14))
                plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='#334155')
                draw_labels(ax, ds, xr.DataArray(data, coords=ds.coords), check_val=15)
                draw_city_labels(ax, xr.DataArray(data, coords=ds.coords), fmt="{:.0f}", suffix=" mph")
                plt.savefig(f"{folder}/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='white'); plt.close()

        # 4. TOTAL PRECIP
        ds_tp = robust_get_data(H, [":APCP:surface"])
        if ds_tp:
            data = gaussian_filter(ds_tp.tp * 0.03937, sigma=SMOOTH_SIGMA)
            fig, ax = setup_plot(H, fxx, "Total Precipitation (in)")
            im = ax.pcolormesh(ds_tp.longitude, ds_tp.latitude, np.where(data > 0.01, data, np.nan),
                               transform=ccrs.PlateCarree(), cmap=ListedColormap(PRECIP_COLORS), norm=BoundaryNorm(PRECIP_LEVELS, len(PRECIP_COLORS)))
            plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='#334155')
            draw_labels(ax, ds_tp, xr.DataArray(data, coords=ds_tp.coords), fmt="{:.2f}", check_val=0.1)
            draw_city_labels(ax, xr.DataArray(data, coords=ds_tp.coords), fmt="{:.2f}", suffix='"')
            plt.savefig(f"frames_total_precip/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='white'); plt.close()

        # 5. PRECIP TYPE (Intensity Based — 5 levels per type)
        ds_pt = robust_get_data(H, [":(REFC|CRAIN|CSNOW|CFRZR|CICEP):"])
        if ds_pt:
            refc = gaussian_filter(ds_pt.refc.values, sigma=SMOOTH_SIGMA)
            # 5 intensity bins: 0-10, 10-20, 20-30, 30-40, 40+ dBZ → indices 0-4
            intensity = np.clip((refc / 10).astype(int), 0, 4)
            final_map = np.zeros_like(refc)
            for i, vname in enumerate(['crain', 'cfrzr', 'icep', 'csnow'], 0):
                if vname in ds_pt: final_map = np.where(ds_pt[vname] == 1, (i*5)+1 + intensity, final_map)
            fig, ax = setup_plot(H, fxx, "Precip Type & Intensity")
            ax.pcolormesh(ds_pt.longitude, ds_pt.latitude, np.ma.masked_where(final_map==0, final_map),
                          transform=ccrs.PlateCarree(), cmap=ListedColormap(PTYPE_COLORS), norm=BoundaryNorm(PTYPE_LEVELS, 21))
            # Custom legend: 4 rows (type) × 5 columns (intensity)
            _types       = [('Rain', 1), ('Frzg Rain', 6), ('Ice Pellets', 11), ('Snow', 16)]
            _intensities = ['V.Light', 'Light', 'Moderate', 'Heavy', 'Intense']
            legend_handles = [
                mpatches.Patch(facecolor=PTYPE_COLORS[base + k], edgecolor='#64748b',
                               linewidth=0.5, label=f'{tname} ({iname})')
                for tname, base in _types
                for k, iname in enumerate(_intensities)
            ]
            ax.legend(handles=legend_handles, ncol=5, loc='lower center',
                      bbox_to_anchor=(0.5, -0.17), fontsize=8, frameon=True,
                      framealpha=0.95, edgecolor='#334155', facecolor='white',
                      columnspacing=0.9, handlelength=1.5, handletextpad=0.4)
            refc_da = xr.DataArray(refc, coords={'latitude': ds_pt.latitude, 'longitude': ds_pt.longitude})
            draw_city_labels(ax, refc_da, fmt="{:.0f}", suffix=" dBZ")
            plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='white'); plt.close()

        # 6. SNOWFALL
        ds_s = robust_get_data(H, [":ASNOW:surface", ":WEASD:surface", ":SNOD:surface"])
        if ds_s:
            data = gaussian_filter(ds_s[list(ds_s.data_vars)[0]] * 39.37, sigma=SMOOTH_SIGMA)
            fig, ax = setup_plot(H, fxx, "Total Snowfall (in)")
            im = ax.pcolormesh(ds_s.longitude, ds_s.latitude, np.where(data > 0.1, data, np.nan),
                               transform=ccrs.PlateCarree(), cmap=ListedColormap(SNOW_COLORS), norm=BoundaryNorm(SNOW_LEVELS, 13))
            plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, aspect=50, shrink=0.8).ax.tick_params(colors='#334155')
            draw_labels(ax, ds_s, xr.DataArray(data, coords=ds_s.coords), fmt="{:.1f}", check_val=0.1)
            draw_city_labels(ax, xr.DataArray(data, coords=ds_s.coords), fmt="{:.1f}", suffix='"')
            plt.savefig(f"frames_snow/f{fxx:02d}.png", dpi=100, bbox_inches='tight', facecolor='white'); plt.close()

        print(f"Finished Frame f{fxx}")
    except Exception as e:
        print(f"Error on f{fxx}: {e}")
