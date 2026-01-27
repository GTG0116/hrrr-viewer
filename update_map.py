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

# Suppress warnings to keep logs clean
warnings.filterwarnings("ignore")

# --- CONFIGURATION ---
# Define the HRRR projection manually as a fallback (Bulletproof fix)
HRRR_CRS = ccrs.LambertConformal(central_longitude=-97.5, 
                                 central_latitude=38.5, 
                                 standard_parallels=(38.5, 38.5))

def add_watermark(ax):
    """Adds the Ephrata Weather logo to the bottom right."""
    ax_inset = ax.inset_axes([0.83, 0.05, 0.14, 0.12])
    ax_inset.set_facecolor('#4da3ff') 
    ax_inset.set_xticks([]); ax_inset.set_yticks([])
    for spine in ax_inset.spines.values():
        spine.set_edgecolor('white'); spine.set_linewidth(2)
    
    ax_inset.text(0.5, 0.75, 'Ephrata', transform=ax_inset.transAxes, 
                  ha='center', va='center', color='white', 
                  fontsize=10, fontweight='bold', family='sans-serif')
    ax_inset.text(0.5, 0.5, '☁️⚡', transform=ax_inset.transAxes, 
                  ha='center', va='center', fontsize=14)
    ax_inset.text(0.5, 0.25, 'Weather', transform=ax_inset.transAxes, 
                  ha='center', va='center', color='white', 
                  fontsize=10, fontweight='bold', family='sans-serif')

def get_precip_cmap():
    """Custom intensity colormap for precip types."""
    colors = ['#ffffff00'] # Transparent
    # Rain (Blue -> Green -> Yellow -> Orange -> Red)
    colors.extend(['#3498db', '#2ecc71', '#f1c40f', '#e67e22', '#e74c3c'])
    # Snow (Dark Blue -> Light Blue)
    colors.extend(['#08306b', '#08519c', '#2171b5', '#6baed6', '#deebf7'])
    # Sleet (Dark Purple -> Light Purple)
    colors.extend(['#3f007d', '#54278f', '#6a51a3', '#9e9ac8', '#efedf5'])
    # Freezing Rain (Dark Pink -> Light Pink)
    colors.extend(['#7a0177', '#ae017e', '#dd3497', '#f768a1', '#fbb4b9'])
    return ListedColormap(colors)

# 1. Find the latest valid run
H_init = None
print("Searching for latest HRRR run...")
for offset in range(24):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        # Check f01 specifically to ensure data exists beyond init
        test_H = Herbie(search_time.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=1)
        if test_H.inventory() is not None:
            H_init = test_H
            print(f"Found valid run: {search_time.strftime('%Y-%m-%d %H:00')}")
            break
    except: continue

if not H_init: 
    print("CRITICAL: No HRRR data found in the last 24 hours.")
    exit(1)

os.makedirs("frames_precip", exist_ok=True)

# 2. Main Generation Loop
for fxx in range(1, 19):
    try:
        # Re-initialize Herbie for the specific forecast hour
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=fxx)
        
        # Download Data (Simulated Reflectivity + P-Types)
        # We use a broad search to ensure we get the variables
        ds = H.xarray(":(REFC|CRAIN|CSNOW|CFRZR|CICEP):")
        if isinstance(ds, list): ds = xr.merge(ds)
        
        # --- FIX: ROBUST CRS DETECTION ---
        # Try to get projection from data, otherwise use the hardcoded HRRR_CRS
        try:
            map_crs = ds.herbie.crs
        except:
            map_crs = HRRR_CRS

        # valid_date handling
        try:
            utc_v = pd.to_datetime(ds.valid_time.values).replace(tzinfo=pytz.UTC)
        except:
            # Fallback if xarray time parsing is weird
            utc_v = H.valid_date.replace(tzinfo=pytz.UTC)
            
        et_v = utc_v.astimezone(pytz.timezone('US/Eastern'))
        timestamp = f"Valid: {et_v.strftime('%m/%d %I:%M %p ET')}"

        # --- DATA PROCESSING ---
        refc = ds.refc.values
        # 5 Intensity bins (0-10, 10-20, etc.)
        intensity = np.clip((refc / 10).astype(int), 0, 4)
        
        # Create Composite Mask
        final_map = np.zeros_like(refc)
        if 'crain' in ds: final_map = np.where(ds.crain == 1, 1 + intensity, final_map)
        if 'csnow' in ds: final_map = np.where(ds.csnow == 1, 6 + intensity, final_map)
        if 'cicep' in ds: final_map = np.where(ds.cicep == 1, 11 + intensity, final_map)
        if 'cfrzr' in ds: final_map = np.where(ds.cfrzr == 1, 16 + intensity, final_map)

        # --- PLOTTING ---
        fig = plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=map_crs)
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.5)
        
        cmap = get_precip_cmap()
        norm = BoundaryNorm(np.arange(0, 22), cmap.N)
        
        ax.pcolormesh(ds.longitude, ds.latitude, np.ma.masked_where(final_map == 0, final_map),
                      transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)

        # Legend
        legend_y = 0.04
        ax.text(0.02, legend_y, "■ RAIN", transform=ax.transAxes, color='#e74c3c', fontweight='bold', fontsize=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        ax.text(0.12, legend_y, "■ SNOW", transform=ax.transAxes, color='#08519c', fontweight='bold', fontsize=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        ax.text(0.22, legend_y, "■ SLEET", transform=ax.transAxes, color='#54278f', fontweight='bold', fontsize=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        ax.text(0.32, legend_y, "■ ICE", transform=ax.transAxes, color='#ae017e', fontweight='bold', fontsize=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

        plt.title(f"EPHRATA WEATHER | Precip Type & Intensity\n{timestamp}", loc='left', fontweight='bold')
        
        add_watermark(ax)
        
        # Save
        plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()
        print(f"Generated frame f{fxx}")

    except Exception as e:
        print(f"Skipping f{fxx} due to error: {e}")
        # We continue the loop so one bad frame doesn't kill the whole viewer
        continue
