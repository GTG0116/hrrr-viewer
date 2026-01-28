import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from herbie import Herbie
import numpy as np
import xarray as xr
from matplotlib.colors import ListedColormap, BoundaryNorm
import os

# --- CONFIGURATION & DIRECTORIES ---
folders = ["frames_temp", "frames_chill", "frames_wind", "frames_gust", "frames_total_precip"]
for f in folders: os.makedirs(f, exist_ok=True)

EXTENTS = [-80.5, -71.5, 38.5, 43.5] 
MAP_CRS = ccrs.LambertConformal(central_longitude=-76.0, central_latitude=41.0)

# --- BRANDED COLOR PALETTES ---

# Wind & Gust Palette (Based on image_ec51d1.png)
WIND_COLORS = [
    '#1a468a', '#2c6eb5', '#4294c7', '#6baed6', '#9ecae1', # 0-40
    '#c6dbef', '#e5f5e0', '#a1d99b', '#74c476', '#31a354', # 40-80
    '#ff9933', '#ff3300', '#990000', '#660000', '#4d004d', # 80-120
    '#ff99cc', '#ff66cc', '#ff33cc', '#ff00cc', '#ff66ff'  # 120-170+
]
WIND_LEVELS = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 200]

# Temp & Chill Palette (Based on image_ec4e17.png)
TEMP_COLORS = [
    '#990033', '#cc0066', '#ff0099', '#ff66cc', '#ffccff', # -60 to -30
    '#e6e6fa', '#ccccff', '#9999ff', '#6666ff', '#0000cc', # -30 to 0
    '#003399', '#0066cc', '#3399ff', '#66ccff', '#99ffff', # 0 to 30
    '#009999', '#006633', '#009933', '#33cc33', '#99ff99', # 30 to 60
    '#ccff99', '#ffff99', '#ffff00', '#ffcc00', '#ff9900', # 60 to 90
    '#ff6600', '#ff3300', '#cc0000', '#990000', '#4d0000'  # 90 to 120+
]
TEMP_LEVELS = list(range(-60, 131, 6)) # 30 colors mapped to 6-degree steps

# Total Precip Palette (NO WHITE - Gap filled with light lime)
PRECIP_COLORS = [
    '#ccff99', '#99ff33', '#33cc33', '#009933', '#006600', # 0.01 - 0.75
    '#004d66', '#0099cc', '#33ccff', '#9999ff', '#9933ff', # 0.75 - 2.0
    '#cc33ff', '#990000', '#cc0000', '#ff3300', '#ff9900'  # 2.0 - 20.0
]
PRECIP_LEVELS = [0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 5.0, 10.0, 15.0, 20.0]

# --- CORE LOGIC ---

def setup_map(title):
    fig = plt.figure(figsize=(12, 10), facecolor='#020617')
    ax = plt.axes(projection=MAP_CRS)
    ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='black', linewidth=1.5, zorder=3)
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='black', linewidth=1.5, zorder=3)
    plt.title(title, loc='left', fontweight='bold', fontsize=14, color='white', pad=10)
    return fig, ax

# 1. Initialize Herbie to find the latest model run
try:
    # Attempt to find the most recent available run
    H_init = Herbie(model='hrrr', product='sfc', priority=['aws', 'nomads'])
    print(f"Using Model Init: {H_init.date}")
except Exception as e:
    print(f"Critical Error: Could not initialize Herbie. {e}")
    exit()

# 2. Process Frames 1-18
for fxx in range(1, 19):
    try:
        # Redefine H for each lead time based on the initialized run
        H = Herbie(H_init.date, model='hrrr', product='sfc', fxx=fxx)
        valid_t = H.valid_date.strftime('%m/%d %I:%M %p')
        
        # --- TEMP & CHILL ---
        for param, folder, name in [("TMP:2 m", "frames_temp", "Temperature"), (":WCHILL:", "frames_chill", "Wind Chill")]:
            ds = H.xarray(param)
            data = (ds[list(ds.data_vars)[0]] - 273.15) * 9/5 + 32
            fig, ax = setup_map(f"HRRR {name} (°F) | {valid_t}")
            ax.pcolormesh(ds.longitude, ds.latitude, data, transform=ccrs.PlateCarree(), 
                          cmap=ListedColormap(TEMP_COLORS), norm=BoundaryNorm(TEMP_LEVELS, len(TEMP_COLORS)))
            plt.savefig(f"{folder}/f{fxx:02d}.png", dpi=90, bbox_inches='tight', facecolor='#020617')
            plt.close()

        # --- WIND & GUST ---
        # Wind requires U and V components; Gust is a single surface variable
        for mode in ['wind', 'gust']:
            if mode == 'wind':
                ds = H.xarray(":(UGRD|VGRD):10 m")
                speed = np.sqrt(ds.u10**2 + ds.v10**2) * 2.237
                name, folder = "Wind Speed", "frames_wind"
            else:
                ds = H.xarray(":GUST:surface")
                speed = ds.gust * 2.237
                name, folder = "Wind Gust", "frames_gust"
            
            fig, ax = setup_map(f"HRRR {name} (mph) | {valid_t}")
            ax.pcolormesh(ds.longitude, ds.latitude, speed, transform=ccrs.PlateCarree(), 
                          cmap=ListedColormap(WIND_COLORS), norm=BoundaryNorm(WIND_LEVELS, len(WIND_COLORS)))
            plt.savefig(f"{folder}/f{fxx:02d}.png", dpi=90, bbox_inches='tight', facecolor='#020617')
            plt.close()

        # --- TOTAL PRECIP ---
        ds = H.xarray(":APCP:surface")
        precip = ds.tp * 0.03937
        fig, ax = setup_map(f"HRRR Total Precip (in) | {valid_t}")
        # Mask out 0 values so background shows through (Clear), then use colors starting at 0.01
        precip_masked = precip.where(precip >= 0.01)
        ax.pcolormesh(ds.longitude, ds.latitude, precip_masked, transform=ccrs.PlateCarree(), 
                      cmap=ListedColormap(PRECIP_COLORS), norm=BoundaryNorm(PRECIP_LEVELS, len(PRECIP_COLORS)))
        plt.savefig(f"frames_total_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight', facecolor='#020617')
        plt.close()

        print(f"Completed Frame f{fxx:02d}")

    except Exception as e:
        print(f"Error processing frame f{fxx:02d}: {e}")
        continue
