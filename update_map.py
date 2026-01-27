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

# --- CUSTOM COLOR SCHEME CONFIG ---
# Intensity Levels (Reflectivity dBZ)
# Rain: Blue (Low) -> Red (High)
# Snow: Dark Blue -> Light Blue
# Sleet: Dark Purple -> Light Purple
# Freezing Rain: Dark Pink -> Light Pink

def get_precip_cmap():
    # 0: None, 1-5: Rain, 6-10: Snow, 11-15: Sleet, 16-20: FRZN Rain
    colors = ['#ffffff00'] # Transparent for no precip
    
    # Rain (Classic Radar: Blue -> Green -> Yellow -> Orange -> Red)
    colors.extend(['#3498db', '#2ecc71', '#f1c40f', '#e67e22', '#e74c3c'])
    
    # Snow (Blue -> Lighter Blue)
    colors.extend(['#0d47a1', '#1976d2', '#42a5f5', '#90caf9', '#e3f2fd'])
    
    # Sleet (Purple -> Lighter Purple)
    colors.extend(['#4a148c', '#7b1fa2', '#9c27b0', '#ba68c8', '#e1bee7'])
    
    # Freezing Rain (Pink -> Lighter Pink)
    colors.extend(['#880e4f', '#c2185b', '#e91e63', '#f06292', '#f8bbd0'])
    
    return ListedColormap(colors)

# 1. Find Data (Starts at f01 to avoid observation/f00 issues)
H_init = None
for offset in range(24):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        # Search specifically for a run that has at least f01 available
        test_H = Herbie(search_time.strftime('%Y-%m-%d %H:00'), model='hrrr', product='sfc', fxx=1)
        if test_H.inventory() is not None:
            H_init = test_H
            break
    except: continue

if not H_init: raise Exception("Data unavailable.")

os.makedirs("frames_precip", exist_ok=True)

# 2. Rendering Loop
for fxx in range(1, 19): # Skip f00, run f01 through f18
    try:
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=fxx)
        
        # Get Data: Reflectivity for intensity + P-Types
        ds = H.xarray(":(REFC|CRAIN|CSNOW|CFRZR|CICEP):")
        if isinstance(ds, list): ds = xr.merge(ds)
        
        refc = ds.refc.values
        
        # Mapping logic:
        # We divide REFC (0-60 dBZ) into 5 intensity bins per type
        intensity = np.clip((refc / 12).astype(int), 0, 4)
        
        final_map = np.zeros_like(refc)
        if 'crain' in ds: final_map = np.where(ds.crain == 1, 1 + intensity, final_map)
        if 'csnow' in ds: final_map = np.where(ds.csnow == 1, 6 + intensity, final_map)
        if 'cicep' in ds: final_map = np.where(ds.cicep == 1, 11 + intensity, final_map) # Sleet
        if 'cfrzr' in ds: final_map = np.where(ds.cfrzr == 1, 16 + intensity, final_map) # FZRA

        fig = plt.figure(figsize=(12, 8), facecolor='white')
        ax = plt.axes(projection=H.crs)
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        
        cmap = get_precip_cmap()
        norm = BoundaryNorm(np.arange(0, 22), cmap.N)
        
        im = ax.pcolormesh(ds.longitude, ds.latitude, final_map, 
                           transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)

        # --- NEW VISUAL LEGEND ---
        # Draw color boxes at the bottom
        leg_y = 0.05
        ax.text(0.02, leg_y, "RAIN", transform=ax.transAxes, color='#e74c3c', fontweight='bold', fontsize=10)
        ax.text(0.12, leg_y, "SNOW", transform=ax.transAxes, color='#0d47a1', fontweight='bold', fontsize=10)
        ax.text(0.22, leg_y, "SLEET", transform=ax.transAxes, color='#4a148c', fontweight='bold', fontsize=10)
        ax.text(0.32, leg_y, "FZ RAIN", transform=ax.transAxes, color='#880e4f', fontweight='bold', fontsize=10)
        ax.text(0.45, leg_y, "Intensity: Darker/Redder = Heavier", transform=ax.transAxes, color='gray', fontsize=8, style='italic')

        plt.title(f"EPHRATA WEATHER | HRRR Precipitation Depiction\nValid: {H.valid_date.strftime('%m/%d %H:%M')}Z", loc='left')
        
        plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
        plt.close()
        print(f"Rendered f{fxx}")

    except Exception as e: print(f"Error f{fxx}: {e}")
