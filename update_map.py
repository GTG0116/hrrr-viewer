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
import warnings

warnings.filterwarnings("ignore")

def add_watermark(ax):
    """Replicates the Ephrata Weather logo in the bottom right corner."""
    # Positioning the logo box in the bottom right
    ax_inset = ax.inset_axes([0.84, 0.05, 0.13, 0.12])
    ax_inset.set_facecolor('#4da3ff') # Light blue from your logo
    ax_inset.set_xticks([]); ax_inset.set_yticks([])
    
    # Border styling
    for spine in ax_inset.spines.values():
        spine.set_edgecolor('white')
        spine.set_linewidth(1.5)
    
    ax_inset.text(0.5, 0.75, 'Ephrata', transform=ax_inset.transAxes, 
                  ha='center', va='center', color='white', 
                  fontsize=9, fontweight='bold')
    ax_inset.text(0.5, 0.5, '☁️', transform=ax_inset.transAxes, 
                  ha='center', va='center', fontsize=12)
    ax_inset.text(0.5, 0.25, 'Weather', transform=ax_inset.transAxes, 
                  ha='center', va='center', color='white', 
                  fontsize=9, fontweight='bold')

# 1. Find the latest data run
H_init = None
for offset in range(24):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        date_str = search_time.strftime('%Y-%m-%d %H:00')
        test_H = Herbie(date_str, model='hrrr', product='sfc', fxx=0, priority=['aws'])
        if test_H.inventory() is not None:
            H_init = test_H
            break
    except: continue

if not H_init: raise Exception("No HRRR data found.")

os.makedirs("frames_temp", exist_ok=True)
os.makedirs("frames_precip", exist_ok=True)

# 2. Render Frames
for fxx in range(19):
    try:
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), model='hrrr', product='sfc', fxx=fxx)
        utc_v = H.valid_date.replace(tzinfo=pytz.UTC)
        et_v = utc_v.astimezone(pytz.timezone('US/Eastern'))
        timestamp = f"Valid: {utc_v.strftime('%m/%d %H:%M')}Z | {et_v.strftime('%I:%M %p ET')}"

        # Get data once to establish grid
        ds_t = H.xarray("TMP:2 m")
        if isinstance(ds_t, list): ds_t = ds_t[0]
        temp_f = (ds_t.t2m - 273.15) * 9/5 + 32
        
        # --- MAP 1: TEMPERATURE ---
        fig = plt.figure(figsize=(12, 7))
        ax = plt.axes(projection=H.crs)
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.5)
        
        im = ax.pcolormesh(ds_t.longitude, ds_t.latitude, temp_f, transform=ccrs.PlateCarree(), 
                           cmap='jet', vmin=-10, vmax=110)
        plt.colorbar(im, label="Temperature (°F)", orientation='horizontal', pad=0.05)
        plt.title(f"HRRR 2m Temperature\n{timestamp}", loc='left', fontweight='bold')
        
        add_watermark(ax)
        plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()

        # --- MAP 2: PRECIP TYPE ---
        ds_p = H.xarray(":(CRAIN|CSNOW|CFRZR|CICEP):")
        precip_mask = np.zeros_like(temp_f)

        def get_mask(d):
            m = np.zeros_like(temp_f)
            if 'crain' in d: m = np.where(d.crain == 1, 1, m)
            if 'csnow' in d: m = np.where(d.csnow == 1, 2, m)
            if 'cfrzr' in d: m = np.where(d.cfrzr == 1, 3, m)
            if 'cicep' in d: m = np.where(d.cicep == 1, 4, m)
            return m

        if isinstance(ds_p, list):
            for sub in ds_p: precip_mask = np.maximum(precip_mask, get_mask(sub))
        else: precip_mask = get_mask(ds_p)

        fig = plt.figure(figsize=(12, 7))
        ax = plt.axes(projection=H.crs)
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        
        # Precip Colors: None, Rain(Grn), Snow(Blu), Ice(Red), Mix(Org)
        p_cmap = ListedColormap(['none', '#2ecc71', '#3498db', '#e74c3c', '#e67e22'])
        ax.pcolormesh(ds_t.longitude, ds_t.latitude, np.ma.masked_where(precip_mask == 0, precip_mask),
                      transform=ccrs.PlateCarree(), cmap=p_cmap, vmin=0, vmax=4)
        
        # Manual Legend for Precip
        plt.text(0.02, 0.05, "■ Rain ■ Snow ■ Ice ■ Mix", transform=ax.transAxes, 
                 fontsize=10, fontweight='bold', color='black', bbox=dict(facecolor='white', alpha=0.7))
        
        plt.title(f"HRRR Precipitation Type\n{timestamp}", loc='left', fontweight='bold')
        add_watermark(ax)
        plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()

        print(f"Generated f{fxx}")
    except Exception as e: print(f"Error at f{fxx}: {e}")
