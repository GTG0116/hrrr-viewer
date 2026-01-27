import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from herbie import Herbie
from datetime import datetime, timedelta
import pytz
import os

# 1. Find the most recent available HRRR run
H = None
for offset in range(4):
    try:
        search_time = datetime.utcnow() - timedelta(hours=offset)
        date_str = search_time.strftime('%Y-%m-%d %H:00')
        test_H = Herbie(date_str, model='hrrr', product='sfc', fxx=0, priority=['aws'])
        if test_H.inventory() is not None:
            H = test_H
            break
    except:
        continue

if not H:
    raise Exception("Could not find recent HRRR data.")

# 2. Determine Forecast Range (48h for big runs, 18h for others)
init_hour = int(H.date.strftime('%H'))
max_fxx = 48 if init_hour in [0, 6, 12, 18] else 18
print(f"Init: {init_hour}z | Max Forecast Available: {max_fxx}h")

# For this viewer, we will grab the Analysis (fxx=0). 
# To show a forecast frame, change fxx to a different number.
ds = H.xarray("TMP:2 m")
temp_f = (ds.t2m - 273.15) * 9/5 + 32

# 3. Time Formatting
utc_time = H.date.replace(tzinfo=pytz.UTC)
eastern_tz = pytz.timezone('US/Eastern')
et_time = utc_time.astimezone(eastern_tz)

init_str = f"Init: {utc_time.strftime('%Y-%m-%d %H:%M')} UTC"
valid_utc = f"Valid: {utc_time.strftime('%Y-%m-%d %H:%M')} UTC"
valid_et = f"({et_time.strftime('%I:%M %p %Z')})"

# 4. Create Plot
fig = plt.figure(figsize=(12, 8))
ax = plt.axes(projection=ds.herbie.crs)
ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='black', linewidth=0.5)

# Plot Data
plot = ax.pcolormesh(ds.longitude, ds.latitude, temp_f, 
                     transform=ccrs.PlateCarree(), cmap='jet', vmin=-20, vmax=110)

plt.colorbar(plot, label="Temperature (°F)", orientation='horizontal', pad=0.05, aspect=50)

# Add Professional Header
plt.title(f"HRRR 2m Temperature\n{init_str} | {valid_utc} {valid_et}", 
          loc='left', fontweight='bold', fontsize=10)

plt.savefig("latest_hrrr.png", dpi=150, bbox_inches='tight')
print("Map updated with time conversions.")
