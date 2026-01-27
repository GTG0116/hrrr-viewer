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

def get_precip_cmap():
    """
    Creates a custom colormap:
    0: Transparent
    1-5: Rain (Blue -> Green -> Yellow -> Orange -> Red)
    6-10: Snow (Dark Blue -> Light Blue)
    11-15: Sleet (Dark Purple -> Light Purple)
    16-20: FZ Rain (Dark Pink -> Light Pink)
    """
    colors = ['#ffffff00'] # 0: None
    
    # Rain (Classic Radar Colors)
    colors.extend(['#3498db', '#2ecc71', '#f1c40f', '#e67e22', '#e74c3c'])
    
    # Snow (Blue Shades - gets lighter)
    colors.extend(['#08306b', '#08519c', '#2171b5', '#6baed6', '#deebf7'])
    
    # Sleet (Purple Shades - gets lighter)
    colors.extend(['#3f007d', '#54278f', '#6a51a3', '#9e9ac8', '#efedf5'])
    
    # Freezing Rain (Pink/Magenta Shades - gets lighter)
    colors.extend(['#7a0177', '#ae017e', '#dd3497', '#f768a1', '#fbb4b9'])
    
    return ListedColormap(colors)

# 1. Find Data (Starts searching at f01 to avoid observation issues)
H_init = None
for offset in range(24):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        test_H = Herbie(search_time.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=1)
        if test_H.inventory() is not None:
            H_init = test_H
            break
    except: continue

if not H_init: raise Exception("HRRR Data unavailable.")

os.makedirs("frames_precip", exist_ok=True)

# 2. Rendering Loop (f01 to f18)
for fxx in range(1, 19):
    try:
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=fxx)
        
        # Load Reflectivity and Precip Types
        ds = H.xarray(":(REFC|CRAIN|CSNOW|CFRZR|CICEP):")
        if isinstance(ds, list): ds = xr.merge(ds)
        
        refc = ds.refc.values
        # Intensity bins: 0-10, 10-20, 20-30, 30-40, 40+ dBZ
        intensity = np.clip((refc / 10).astype(int), 0, 4)
        
        # Combine masks with intensity
        # Rain=1, Snow=6, Sleet=11, FZRA=16 (+ intensity 0-4)
        final_map = np.zeros_like(refc)
        if 'crain' in ds: final_map = np.where(ds.crain == 1, 1 + intensity, final_map)
        if 'csnow' in ds: final_map = np.where(ds.csnow == 1, 6 + intensity, final_map)
        if 'cicep' in ds: final_map = np.where(ds.cicep == 1, 11 + intensity, final_map)
        if 'cfrzr' in ds: final_map = np.where(ds.cfrzr == 1, 16 + intensity, final_map)

        fig = plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=H.crs)
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        
        cmap = get_precip_cmap()
        norm = BoundaryNorm(np.arange(0, 22), cmap.N)
        
        ax.pcolormesh(ds.longitude, ds.latitude, np.ma.masked_where(final_map == 0, final_map),
                      transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)

        # Legend with matching colors
        legend_data = [
            ("RAIN", "#e74c3c", 0.02),
            ("SNOW", "#08519c", 0.12),
            ("SLEET", "#54278f", 0.22),
            ("FZ RAIN", "#ae017e", 0.32)
        ]
        for text, color, x_pos in legend_data:
            ax.text(x_pos, 0.05, f"■ {text}", transform=ax.transAxes, 
                    color=color, fontweight='bold', fontsize=12,
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

        plt.title(f"EPHRATA WEATHER | HRRR Precip Depiction\nValid: {H.valid_date.strftime('%m/%d %H:%M')}Z", loc='left', fontweight='bold')
        
        plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()
        print(f"Success: f{fxx}")
        
    except Exception as e: print(f"Error f{fxx}: {e}")
