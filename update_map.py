import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from herbie import Herbie
from datetime import datetime, timedelta
import os

# 1. Logic to find the most recent HRRR data
try:
    # Try the current hour
    now = datetime.utcnow()
    date_str = now.strftime('%Y-%m-%d %H:00')
    print(f"Attempting to fetch data for: {date_str}")
    H = Herbie(date_str, model='hrrr', product='sfc', fxx=0, priority=['aws'])
except Exception as e:
    # If the current hour isn't on AWS yet, try the previous hour
    print(f"Current hour not found, trying previous hour. Error: {e}")
    yesterday = datetime.utcnow() - timedelta(hours=1)
    date_str = yesterday.strftime('%Y-%m-%d %H:00')
    H = Herbie(date_str, model='hrrr', product='sfc', fxx=0, priority=['aws'])

# 2. Decode the GRIB2 data
ds = H.xarray("TMP:2 m")
temp_f = (ds.t2m - 273.15) * 9/5 + 32

# 3. Create the Visualization
fig = plt.figure(figsize=(12, 7))
ax = plt.axes(projection=ds.herbie.crs)
ax.coastlines(resolution='50m', color='black', linewidth=1)

# Add temperature data
plot = ax.pcolormesh(ds.longitude, ds.latitude, temp_f, 
                     transform=ccrs.PlateCarree(), cmap='RdYlBu_r')

plt.colorbar(plot, label="Temperature (°F)", orientation='vertical', pad=0.02)
plt.title(f"HRRR 2m Temperature - Valid: {H.date}")

# 4. Save the final image
plt.savefig("latest_hrrr.png", dpi=150, bbox_inches='tight')
print("Map generated and saved as latest_hrrr.png")
