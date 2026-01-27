import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from herbie import Herbie
from datetime import datetime, timedelta
import pytz
import os
import numpy as np
import xarray as xr
from matplotlib.colors import ListedColormap
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import warnings

warnings.filterwarnings("ignore")

def add_watermark(ax):
    """Adds a replica of the Ephrata Weather logo in the corner."""
    # Create light blue box for the logo
    ax_inset = ax.inset_axes([0.83, 0.05, 0.14, 0.12])
    ax_inset.set_facecolor('#4da3ff') 
    ax_inset.set_xticks([]); ax_inset.set_yticks([])
    for spine in ax_inset.spines.values():
        spine.set_edgecolor('white')
        spine.set_linewidth(2)
    
    ax_inset.text(0.5, 0.75, 'Ephrata', transform=ax_inset.transAxes, 
                  ha='center', va='center', color='white', 
                  fontsize=10, fontweight='bold', family='sans-serif')
    ax_inset.text(0.5, 0.5, '☁️⚡', transform=ax_inset.transAxes, 
                  ha='center', va='center', fontsize=14)
    ax_inset.text(0.5, 0.25, 'Weather', transform=ax_inset.transAxes, 
                  ha='center', va='center', color='white', 
                  fontsize=10, fontweight='bold', family='sans-serif')

# 1. Find the latest available data
H_init = None
for offset in range(24):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        date_str = search_time.strftime('%Y-%m-%d %H:00')
        test_H = Herbie(date_str, model='hrrr', product='sfc', fxx=0, priority=['aws'])
        if test_H.inventory() is not None:
            H_init = test_H
            print(f"Success! Found run: {date_str}")
            break
    except: continue

if not H_init: raise Exception("No HRRR data available.")

os.makedirs("frames_temp", exist_ok=True)
os.makedirs("frames_precip", exist_ok=True)

# 2. Render Frames
for fxx in range(19): # 0 to 18 hours
    try:
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), model='hrrr', product='sfc', fxx=fxx)
        utc_v = H.valid_date.replace(tzinfo=pytz.UTC)
        et_v = utc_v.astimezone(pytz.timezone('US/Eastern'))
        timestamp = f"Valid: {utc_v.strftime('%m/%d %H:%M')}Z | {et_v.strftime('%I:%M %p ET')}"

        # Fetch Temperature and Get Projection from Data
        ds_t = H.xarray("TMP:2 m")
        if isinstance(ds_t, list): ds_t = ds_t[0]
        
        # FIX: Robustly get CRS from the dataset metadata
        try:
            map_crs = ds_t.herbie.crs
        except:
            map_crs = ccrs.LambertConformal(central_longitude=-97.5, central_latitude=38.5)

        temp_f = (ds_t.t2m - 273.15) * 9/5 + 32

        # --- TEMPERATURE MAP ---
        fig = plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=map_crs)
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.5)
        
        im = ax.pcolormesh(ds_t.longitude, ds_t.latitude, temp_f, transform=ccrs.PlateCarree(), 
                           cmap='jet', vmin=-10, vmax=110)
        plt.colorbar(im, label="Temperature (°F)", orientation='horizontal', pad=0.05)
        plt.title(f"HRRR 2m Temperature\n{timestamp}", loc='left', fontweight='bold')
        
        add_watermark(ax)
        plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()

        # --- PRECIPITATION TYPE MAP ---
        ds_p = H.xarray(":(CRAIN|CSNOW|CFRZR|CICEP):")
        precip_mask = np.zeros_like(temp_f)
        
        def process_mask(d):
            m = np.zeros_like(temp_f)
            if 'crain' in d: m = np.where(d.crain == 1, 1, m)
            if 'csnow' in d: m = np.where(d.csnow == 1, 2, m)
            if 'cfrzr' in d: m = np.where(d.cfrzr == 1, 3, m)
            if 'cicep' in d: m = np.where(d.cicep == 1, 4, m)
            return m

        if isinstance(ds_p, list):
            for sub in ds_p: precip_mask = np.maximum(precip_mask, process_mask(sub))
        else: precip_mask = process_mask(ds_p)

        fig = plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=map_crs)
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        
        p_cmap = ListedColormap(['none', '#2ecc71', '#3498db', '#e74c3c', '#e67e22'])
        ax.pcolormesh(ds_t.longitude, ds_t.latitude, np.ma.masked_where(precip_mask == 0, precip_mask),
                      transform=ccrs.PlateCarree(), cmap=p_cmap, vmin=0, vmax=4)
        
        # Legend for Precip
        plt.text(0.02, 0.05, "■ Rain ■ Snow ■ Ice ■ Mix", transform=ax.transAxes, 
                 fontsize=11, fontweight='bold', bbox=dict(facecolor='white', alpha=0.8))
        
        plt.title(f"HRRR Precipitation Type\n{timestamp}", loc='left', fontweight='bold')
        add_watermark(ax)
        plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()
        
        print(f"Finished hour f{fxx}")
    except Exception as e: print(f"Error at f{fxx}: {e}")
