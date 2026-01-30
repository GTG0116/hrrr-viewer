import os
import gzip
import shutil
import requests
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from datetime import datetime, timedelta
import numpy as np

# --- CONFIGURATION (Match your HRRR script) ---
EXTENTS = [-80.5, -71.5, 38.5, 43.5] 
MAP_CRS = ccrs.LambertConformal(central_longitude=-76.0, central_latitude=41.0)
NUM_SCANS = 10  # Last 10 scans

# Folders
os.makedirs("frames_mrms_refl", exist_ok=True)
os.makedirs("frames_mrms_ptype", exist_ok=True)

def fetch_mrms_latest(product, count=10):
    """Fetches the last X timestamps for a specific MRMS product."""
    base_url = "https://noaa-mrms-pds.s3.amazonaws.com/CONUS"
    now = datetime.utcnow()
    # MRMS updates every 2 mins; start looking 4 mins ago to ensure availability
    current_time = now - timedelta(minutes=4)
    found_files = []
    
    while len(found_files) < count:
        timestamp = current_time.strftime("%Y%m%d-%H%M")
        date_path = current_time.strftime("%Y%m%d")
        # Note: filenames use SS (seconds), often 00, but can vary. 
        # For simplicity, we search for HHMM00.
        filename = f"MRMS_{product}_00.00_{timestamp}00.grib2.gz"
        url = f"{base_url}/{product}/{date_path}/{filename}"
        
        r = requests.head(url)
        if r.status_code == 200:
            found_files.append((url, filename, current_time))
        
        current_time -= timedelta(minutes=2)
        if (now - current_time).total_seconds() > 3600: # Stop after 1 hour of missing data
            break
            
    return sorted(found_files) # Oldest first

def process_mrms():
    # Fetch lists
    refl_files = fetch_mrms_latest("LayerCompositeReflectivity_Low", NUM_SCANS)
    ptype_files = fetch_mrms_latest("PrecipFlag", NUM_SCANS)

    for i, (url, fname, ts) in enumerate(refl_files):
        fxx = i + 1
        print(f"Processing MRMS Frame {fxx}/10: {ts.strftime('%H:%M')}Z")
        
        # Download and Decompress
        local_gz = f"temp_{fxx}.grib2.gz"
        local_grib = f"temp_{fxx}.grib2"
        with requests.get(url, stream=True) as r:
            with open(local_gz, 'wb') as f:
                shutil.copyfileobj(r.raw, f)
        
        with gzip.open(local_gz, 'rb') as f_in:
            with open(local_grib, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        # Decode and Plot Reflectivity
        try:
            ds = xr.open_dataset(local_grib, engine='cfgrib')
            refl = ds.unknown # Variable name is often 'unknown' in these GRIBs
            
            fig = plt.figure(figsize=(10, 8), facecolor='black')
            ax = plt.axes(projection=MAP_CRS)
            ax.set_extent(EXTENTS, crs=ccrs.PlateCarree())
            
            # Use standard NWS Reflectivity colors
            levels = np.arange(5, 75, 5)
            im = refl.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), levels=levels, 
                                     cmap='pyart_NWSRef', add_colorbar=False)
            
            # Save
            plt.savefig(f"frames_mrms_refl/f{fxx:02d}.png", dpi=90, bbox_inches='tight')
            plt.close()
        except Exception as e: print(f"Error {fxx}: {e}")
        
        # Cleanup
        if os.path.exists(local_gz): os.remove(local_gz)
        if os.path.exists(local_grib): os.remove(local_grib)

if __name__ == "__main__":
    process_mrms()
