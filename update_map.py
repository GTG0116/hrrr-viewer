import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from herbie import Herbie
from datetime import datetime, timedelta
import os

# 1. Logic to find the most RECENT available HRRR data
H = None
# We will try the current hour, then 1 hour ago, then 2 hours ago
for offset in range(3):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        date_str = search_time.strftime('%Y-%m-%d %H:00')
        print(f"Checking for data at: {date_str}...")
        
        test_H = Herbie(date_str, model='hrrr', product='sfc', fxx=0, priority=['aws'])
        
        # Try to load the index to verify it exists
        if test_H.inventory() is not None:
            H = test_H
            print(f"Success! Using data from: {date_str}")
            break
    except Exception as e:
        print(f"Data not ready for {date_str}, trying previous hour...")

if H is None:
    raise Exception("Could not find any HRRR data from the last 3 hours.")

# 2. Decode the GRIB2 data
ds = H.xarray("TMP:2 m")
temp_f = (ds.t2m - 273.15) * 9/5 + 32

# 3. Create the Visualization
fig = plt.figure(figsize=(12, 7))
ax = plt.axes(projection=ds.herbie.crs)
ax.coastlines(resolution='50m', color='black', linewidth=1)

# This adds state borders to make the map useful
import cartopy.feature as cfeature
ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='gray', linewidth=0.5)

# Add temperature data
plot = ax.pcolormesh(ds.longitude, ds.latitude, temp_f, 
                     transform=ccrs.PlateCarree(), cmap='jet')

plt.colorbar(plot, label="Temperature (°F)", orientation='vertical', pad=0.02)
plt.title(f"HRRR 2m Temperature - Valid: {H.date} UTC")

# 4. Save the final image
plt.savefig("latest_hrrr.png", dpi=150, bbox_inches='tight')
print("Map successfully saved as latest_hrrr.png")
