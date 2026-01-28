import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
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
# Northeast Extents (approximate)
EXTENTS = [-82, -67, 37, 48] # [West, East, South, North]

# Custom Projection (Lambert Conformal focused on NE)
# We define this manually to avoid the "no attribute 'crs'" error
MAP_CRS = ccrs.LambertConformal(central_longitude=-74.5, central_latitude=42.0)

def add_watermark(ax):
    """Adds the Ephrata Weather logo."""
    # Adjusted position for the new Northeast zoom
    ax_inset = ax.inset_axes([0.80, 0.02, 0.18, 0.12])
    ax_inset.set_facecolor('#4da3ff') 
    ax_inset.set_xticks([]); ax_inset.set_yticks([])
    for spine in ax_inset.spines.values():
        spine.set_edgecolor('white'); spine.set_linewidth(1.5)
    
    ax_inset.text(0.5, 0.75, 'Ephrata', transform=ax_inset.transAxes, 
                  ha='center', va='center', color='white', 
                  fontsize=9, fontweight='bold', family='sans-serif')
    ax_inset.text(0.5, 0.5, '☁️⚡', transform=ax_inset.transAxes, 
                  ha='center', va='center', fontsize=12)
    ax_inset.text(0.5, 0.25, 'Weather', transform=ax_inset.transAxes, 
                  ha='center', va='center', color='white', 
                  fontsize=9, fontweight='bold', family='sans-serif')

def setup_map(title):
    """Helper to create a fresh figure with common settings."""
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=MAP_CRS)
    ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='black', linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='black', linewidth=0.7)
    add_watermark(ax)
    plt.title(title, loc='left', fontweight='bold', fontsize=10)
    return fig, ax

def get_precip_cmap():
    """Custom intensity colormap."""
    colors = ['#ffffff00'] 
    colors.extend(['#3498db', '#2ecc71', '#f1c40f', '#e67e22', '#e74c3c']) # Rain
    colors.extend(['#08306b', '#08519c', '#2171b5', '#6baed6', '#deebf7']) # Snow
    colors.extend(['#3f007d', '#54278f', '#6a51a3', '#9e9ac8', '#efedf5']) # Sleet
    colors.extend(['#7a0177', '#ae017e', '#dd3497', '#f768a1', '#fbb4b9']) # Ice
    return ListedColormap(colors)

# 1. Find Data
H_init = None
print("Searching for latest HRRR run...")
for offset in range(24):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        # Check f01 to ensure data exists
        test_H = Herbie(search_time.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=1)
        if test_H.inventory() is not None:
            H_init = test_H
            print(f"Found run: {search_time.strftime('%Y-%m-%d %H:00')}")
            break
    except: continue

if not H_init: exit(1)

# Create directories
folders = ["frames_temp", "frames_precip", "frames_wind", "frames_gust", "frames_snow", "frames_total_precip"]
for folder in folders:
    os.makedirs(folder, exist_ok=True)

# 2. Main Generation Loop
for fxx in range(1, 19):
    try:
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=fxx)
        
        # Get standardized time strings
        utc_v = H.valid_date.replace(tzinfo=pytz.UTC)
        et_v = utc_v.astimezone(pytz.timezone('US/Eastern'))
        t_str = f"Valid: {et_v.strftime('%m/%d %I:%M %p ET')}"

        # ---------------------------------------------------------
        # 1. TEMPERATURE (2m)
        # ---------------------------------------------------------
        try:
            ds_t = H.xarray("TMP:2 m")
            if isinstance(ds_t, list): ds_t = ds_t[0]
            temp_f = (ds_t.t2m - 273.15) * 9/5 + 32
            
            fig, ax = setup_map(f"HRRR 2m Temperature | {t_str}")
            im = ax.pcolormesh(ds_t.longitude, ds_t.latitude, temp_f, transform=ccrs.PlateCarree(), 
                               cmap='jet', vmin=-10, vmax=100)
            plt.colorbar(im, label="Temp (°F)", orientation='horizontal', pad=0.02, shrink=0.8)
            plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Temp Error f{fxx}: {e}")

        # ---------------------------------------------------------
        # 2. PRECIPITATION TYPE & INTENSITY
        # ---------------------------------------------------------
        try:
            ds_p = H.xarray(":(REFC|CRAIN|CSNOW|CFRZR|CICEP):")
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
            norm = BoundaryNorm(np.arange(0, 22), cmap.N)
            ax.pcolormesh(ds_p.longitude, ds_p.latitude, np.ma.masked_where(final_map==0, final_map),
                          transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
            
            # Simple Legend
            ax.text(0.02, 0.02, "■ RAIN  ■ SNOW  ■ SLEET  ■ ICE", transform=ax.transAxes, 
                    fontsize=8, fontweight='bold', bbox=dict(facecolor='white', alpha=0.8))
            
            plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Precip Error f{fxx}: {e}")

        # ---------------------------------------------------------
        # 3. WIND (10m)
        # ---------------------------------------------------------
        try:
            ds_w = H.xarray(":(UGRD|VGRD):10 m")
            if isinstance(ds_w, list): ds_w = xr.merge(ds_w)
            # Calculate magnitude (speed) in knots (approx 1.944 * m/s)
            wind_speed = np.sqrt(ds_w.u10**2 + ds_w.v10**2) * 1.944
            
            fig, ax = setup_map(f"HRRR 10m Wind Speed | {t_str}")
            im = ax.pcolormesh(ds_w.longitude, ds_w.latitude, wind_speed, transform=ccrs.PlateCarree(),
                               cmap='BuPu', vmin=0, vmax=50)
            plt.colorbar(im, label="Wind Speed (knots)", orientation='horizontal', pad=0.02, shrink=0.8)
            plt.savefig(f"frames_wind/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Wind Error f{fxx}: {e}")

        # ---------------------------------------------------------
        # 4. WIND GUSTS (Surface)
        # ---------------------------------------------------------
        try:
            ds_g = H.xarray("GUST:surface")
            if isinstance(ds_g, list): ds_g = ds_g[0]
            gust_kts = ds_g.gust * 1.944
            
            fig, ax = setup_map(f"HRRR Wind Gusts | {t_str}")
            im = ax.pcolormesh(ds_g.longitude, ds_g.latitude, gust_kts, transform=ccrs.PlateCarree(),
                               cmap='Reds', vmin=0, vmax=60)
            plt.colorbar(im, label="Gusts (knots)", orientation='horizontal', pad=0.02, shrink=0.8)
            plt.savefig(f"frames_gust/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Gust Error f{fxx}: {e}")

        # ---------------------------------------------------------
        # 5. TOTAL SNOWFALL (Accumulated)
        # ---------------------------------------------------------
        try:
            ds_s = H.xarray("ASNOW:surface")
            if isinstance(ds_s, list): ds_s = ds_s[0]
            # Convert meters to inches
            snow_in = ds_s.asnow * 39.37
            
            fig, ax = setup_map(f"HRRR Total Snowfall | {t_str}")
            im = ax.pcolormesh(ds_s.longitude, ds_s.latitude, snow_in, transform=ccrs.PlateCarree(),
                               cmap='cool', vmin=0, vmax=12) # 0-12 inches scale
            plt.colorbar(im, label="Accumulated Snow (in)", orientation='horizontal', pad=0.02, shrink=0.8)
            plt.savefig(f"frames_snow/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Snow Error f{fxx}: {e}")

        # ---------------------------------------------------------
        # 6. TOTAL PRECIPITATION (Accumulated)
        # ---------------------------------------------------------
        try:
            ds_tp = H.xarray("APCP:surface")
            if isinstance(ds_tp, list): ds_tp = ds_tp[0]
            # Convert kg/m^2 (mm) to inches
            precip_in = ds_tp.tp * 0.03937
            
            fig, ax = setup_map(f"HRRR Total Precipitation | {t_str}")
            im = ax.pcolormesh(ds_tp.longitude, ds_tp.latitude, precip_in, transform=ccrs.PlateCarree(),
                               cmap='gist_earth_r', vmin=0, vmax=3) # 0-3 inches scale
            plt.colorbar(im, label="Total Precip (in)", orientation='horizontal', pad=0.02, shrink=0.8)
            plt.savefig(f"frames_total_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Total Precip Error f{fxx}: {e}")

        print(f"Completed Frame f{fxx}")

    except Exception as e:
        print(f"Skipping f{fxx} entirely: {e}")
