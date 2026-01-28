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
import warnings
from shapely.geometry import shape

warnings.filterwarnings("ignore")

# --- CONFIGURATION ---
EXTENTS = [-80.5, -71.5, 38.5, 43.5] # Focused Northeast (PA/NY/NJ/MD/DE/CT area)
# We zoom in slightly more than before so the county numbers are readable

MAP_CRS = ccrs.LambertConformal(central_longitude=-76.0, central_latitude=41.0)

# Pre-load County Shapes (Efficiency)
print("Pre-loading US Counties...")
reader = shpreader.Reader(shpreader.natural_earth(resolution='10m', category='cultural', name='admin_2_counties'))
counties = list(reader.records())
print(f"Loaded {len(counties)} counties.")

def get_main_var(ds):
    """Helper to safely get the first data variable from a dataset."""
    if isinstance(ds, list):
        # If list, try to merge, or take the first one that has data
        try: ds = xr.merge(ds, compat='override')
        except: ds = ds[0]
    
    # Return the values of the first variable found
    var_name = list(ds.data_vars)[0]
    return ds, ds[var_name]

def draw_county_labels(ax, ds, data_var, fmt="{:.0f}", color='black', check_val=None):
    """
    Iterates through counties, finds centroid, samples data, and prints value.
    check_val: If provided, only print if value > check_val (good for snow/precip)
    """
    # Filter counties strictly within our view to save time
    view_w, view_e, view_s, view_n = EXTENTS
    
    for county in counties:
        geom = county.geometry
        bounds = geom.bounds # (minx, miny, maxx, maxy)
        
        # Quick bounding box check
        if bounds[2] < view_w or bounds[0] > view_e or bounds[3] < view_s or bounds[1] > view_n:
            continue
            
        centroid = geom.centroid
        lon, lat = centroid.x, centroid.y
        
        # Sample the data at this point (Nearest Neighbor)
        try:
            # We select by lat/lon. Note: HRRR uses grid, so we use the coordinates provided
            # Xarray selection requires projection handling if coords are x/y, 
            # but usually Herbie returns lat/lon coords we can use if we are careful.
            # A safer generic way for unstructured grids is finding the index:
            
            # Simple nearest neighbor using xarray's selection if coordinates allow
            if 'latitude' in ds.coords and 'longitude' in ds.coords:
                # Calculate distance to find nearest point index manually (most robust for curvelinear grids)
                # This is computationally expensive for every county. 
                # OPTIMIZATION: simple box selection first? 
                # Let's use the provided .sel method if lat/lon are dimensions, otherwise fallback.
                val = data_var.sel(latitude=lat, longitude=lon, method='nearest').values
            else:
                # Fallback for projection coordinates (y, x)
                # This is complex without reprojecting. 
                # Simplified approach: skip readout if direct select fails to avoid crashes
                continue
                
            val = float(val)
            
            # Filter zeros for things like snow/precip (don't print "0" everywhere)
            if check_val is not None and val < check_val:
                continue

            # Plot text
            ax.text(lon, lat, fmt.format(val), transform=ccrs.PlateCarree(),
                    ha='center', va='center', fontsize=6, fontweight='bold',
                    color=color, clip_on=True, zorder=10)
        except:
            continue

def add_watermark(ax):
    """Adds the Ephrata Weather logo."""
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
    fig = plt.figure(figsize=(12, 10)) # Larger figure for labels
    ax = plt.axes(projection=MAP_CRS)
    ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
    
    # Add Features
    ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='black', linewidth=1.5, zorder=2)
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='black', linewidth=1, zorder=2)
    
    # Add Counties (Light Grey)
    ax.add_feature(cfeature.ShapelyFeature(shpreader.Reader(shpreader.natural_earth(resolution='10m', category='cultural', name='admin_2_counties')).geometries(), ccrs.PlateCarree()), 
                   facecolor='none', edgecolor='gray', linewidth=0.3, zorder=1)

    add_watermark(ax)
    plt.title(title, loc='left', fontweight='bold', fontsize=12)
    return fig, ax

def get_precip_cmap():
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
        test_H = Herbie(search_time.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=1)
        if test_H.inventory() is not None:
            H_init = test_H
            print(f"Found run: {search_time.strftime('%Y-%m-%d %H:00')}")
            break
    except: continue

if not H_init: exit(1)

folders = ["frames_temp", "frames_precip", "frames_wind", "frames_gust", "frames_snow", "frames_total_precip"]
for folder in folders:
    os.makedirs(folder, exist_ok=True)

# 2. Main Generation Loop
for fxx in range(1, 19):
    try:
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), model='hrrr', fxx=fxx)
        
        utc_v = H.valid_date.replace(tzinfo=pytz.UTC)
        et_v = utc_v.astimezone(pytz.timezone('US/Eastern'))
        t_str = f"Valid: {et_v.strftime('%m/%d %I:%M %p ET')}"

        # ---------------------------------------------------------
        # 1. TEMPERATURE (2m)
        # ---------------------------------------------------------
        try:
            # Use broader search to ensure we get it
            ds_raw = H.xarray("TMP:2 m", verbose=False)
            ds_t, temp_var = get_main_var(ds_raw)
            
            # Convert K to F
            temp_f = (temp_var - 273.15) * 9/5 + 32
            
            fig, ax = setup_map(f"HRRR 2m Temperature (°F) | {t_str}")
            im = ax.pcolormesh(ds_t.longitude, ds_t.latitude, temp_f, transform=ccrs.PlateCarree(), 
                               cmap='jet', vmin=-10, vmax=100)
            
            # Readouts (No threshold, print all)
            # Since temp needs conversion, we pass the converted array
            # Note: The 'draw_county_labels' needs a DataArray. We create a temp one.
            temp_f_da = xr.DataArray(temp_f, coords=ds_t.coords, dims=ds_t.dims)
            draw_county_labels(ax, ds_t, temp_f_da, fmt="{:.0f}", color='black')

            plt.colorbar(im, fraction=0.046, pad=0.04)
            plt.savefig(f"frames_temp/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Temp Error f{fxx}: {e}")

        # ---------------------------------------------------------
        # 2. PRECIPITATION TYPE
        # ---------------------------------------------------------
        try:
            ds_p = H.xarray(":(REFC|CRAIN|CSNOW|CFRZR|CICEP):", verbose=False)
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
            
            # No text readouts for this map (too cluttered with p-type codes)
            
            plt.savefig(f"frames_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Precip Error f{fxx}: {e}")

        # ---------------------------------------------------------
        # 3. WIND SPEED
        # ---------------------------------------------------------
        try:
            ds_w = H.xarray(":(UGRD|VGRD):10 m", verbose=False)
            if isinstance(ds_w, list): ds_w = xr.merge(ds_w)
            
            ws = np.sqrt(ds_w.u10**2 + ds_w.v10**2) * 1.944 # to knots
            
            fig, ax = setup_map(f"HRRR Wind Speed (kts) | {t_str}")
            im = ax.pcolormesh(ds_w.longitude, ds_w.latitude, ws, transform=ccrs.PlateCarree(),
                               cmap='BuPu', vmin=0, vmax=50)
            
            ws_da = xr.DataArray(ws, coords=ds_w.coords, dims=ds_w.dims)
            draw_county_labels(ax, ds_w, ws_da, fmt="{:.0f}", color='black', check_val=5) # Only show > 5 kts

            plt.colorbar(im, fraction=0.046, pad=0.04)
            plt.savefig(f"frames_wind/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Wind Error f{fxx}: {e}")

        # ---------------------------------------------------------
        # 4. WIND GUSTS
        # ---------------------------------------------------------
        try:
            ds_raw = H.xarray(":GUST:surface", verbose=False)
            ds_g, gust_var = get_main_var(ds_raw)
            gust_kts = gust_var * 1.944
            
            fig, ax = setup_map(f"HRRR Wind Gusts (kts) | {t_str}")
            im = ax.pcolormesh(ds_g.longitude, ds_g.latitude, gust_kts, transform=ccrs.PlateCarree(),
                               cmap='Reds', vmin=0, vmax=60)
            
            gust_da = xr.DataArray(gust_kts, coords=ds_g.coords, dims=ds_g.dims)
            draw_county_labels(ax, ds_g, gust_da, fmt="{:.0f}", color='black', check_val=10) # Only show > 10

            plt.colorbar(im, fraction=0.046, pad=0.04)
            plt.savefig(f"frames_gust/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Gust Error f{fxx}: {e}")

        # ---------------------------------------------------------
        # 5. TOTAL SNOWFALL
        # ---------------------------------------------------------
        try:
            # HRRR uses ASNOW for accumulation. Sometimes it's WEASD in other models.
            ds_raw = H.xarray(":ASNOW:surface", verbose=False)
            ds_s, snow_var = get_main_var(ds_raw)
            snow_in = snow_var * 39.37
            
            fig, ax = setup_map(f"HRRR Total Snow (in) | {t_str}")
            im = ax.pcolormesh(ds_s.longitude, ds_s.latitude, snow_in, transform=ccrs.PlateCarree(),
                               cmap='cool', vmin=0, vmax=12)
            
            snow_da = xr.DataArray(snow_in, coords=ds_s.coords, dims=ds_s.dims)
            draw_county_labels(ax, ds_s, snow_da, fmt="{:.1f}", color='black', check_val=0.1) # Only > 0.1 inch

            plt.colorbar(im, fraction=0.046, pad=0.04)
            plt.savefig(f"frames_snow/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Snow Error f{fxx}: {e}")

        # ---------------------------------------------------------
        # 6. TOTAL PRECIP
        # ---------------------------------------------------------
        try:
            ds_raw = H.xarray(":APCP:surface", verbose=False)
            ds_tp, precip_var = get_main_var(ds_raw)
            precip_in = precip_var * 0.03937
            
            fig, ax = setup_map(f"HRRR Total Precip (in) | {t_str}")
            im = ax.pcolormesh(ds_tp.longitude, ds_tp.latitude, precip_in, transform=ccrs.PlateCarree(),
                               cmap='gist_earth_r', vmin=0, vmax=3)
            
            precip_da = xr.DataArray(precip_in, coords=ds_tp.coords, dims=ds_tp.dims)
            draw_county_labels(ax, ds_tp, precip_da, fmt="{:.2f}", color='blue', check_val=0.01)

            plt.colorbar(im, fraction=0.046, pad=0.04)
            plt.savefig(f"frames_total_precip/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Total Precip Error f{fxx}: {e}")

        print(f"Completed Frame f{fxx}")

    except Exception as e:
        print(f"Skipping f{fxx} entirely: {e}")
