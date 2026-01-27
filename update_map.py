import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from herbie import Herbie
from datetime import datetime, timedelta
import pytz
import os

# 1. Find the most recent available HRRR run
H_init = None
for offset in range(4):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        date_str = search_time.strftime('%Y-%m-%d %H:00')
        test_H = Herbie(date_str, model='hrrr', product='sfc', fxx=0, priority=['aws'])
        if test_H.inventory() is not None:
            H_init = test_H
            break
    except:
        continue

if not H_init:
    raise Exception("Could not find recent HRRR data.")

# 2. Determine Forecast Range
init_hour = int(H_init.date.strftime('%H'))
max_fxx = 48 if init_hour in [0, 6, 12, 18] else 18
print(f"Init: {init_hour}z | Generating {max_fxx} forecast frames.")

# Create a folder for frames if it doesn't exist
os.makedirs("frames", exist_ok=True)

# 3. Loop through forecast hours
for fxx in range(max_fxx + 1):
    try:
        # Fetch specific forecast hour
        H = Herbie(H_init.date.strftime('%Y-%m-%d %H:00'), 
                   model='hrrr', product='sfc', fxx=fxx, priority=['aws'])
        
        ds = H.xarray("TMP:2 m")
        temp_f = (ds.t2m - 273.15) * 9/5 + 32

        # Time Math
        utc_valid = H.valid_date.replace(tzinfo=pytz.UTC)
        et_valid = utc_valid.astimezone(pytz.timezone('US/Eastern'))

        # Plotting
        fig = plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=ds.herbie.crs)
        ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='black', linewidth=0.5)
        
        plot = ax.pcolormesh(ds.longitude, ds.latitude, temp_f, 
                             transform=ccrs.PlateCarree(), cmap='jet', vmin=-20, vmax=110)
        
        plt.colorbar(plot, label="Temp (°F)", orientation='horizontal', pad=0.05, aspect=50)

        # Header with Init and Valid times
        title_str = (f"HRRR Init: {H.date.strftime('%m/%d %H:%M')} UTC | "
                     f"Valid: {utc_valid.strftime('%m/%d %H:%M')} UTC\n"
                     f"Forecast Hour: f{fxx:02d} | Eastern: {et_valid.strftime('%I:%M %p %Z')}")
        plt.title(title_str, loc='left', fontweight='bold', fontsize=10)

        # Save individual frame
        plt.savefig(f"frames/frame_{fxx:02d}.png", dpi=100, bbox_inches='tight')
        plt.close() # Close plot to save memory
        print(f"Generated frame {fxx}")

    except Exception as e:
        print(f"Skipping fxx {fxx}: {e}")
