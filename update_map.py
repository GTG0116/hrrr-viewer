import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from herbie import Herbie  # <--- Make sure this stays 'herbie'
import os

# 1. Fetch the data (fxx=0 is the most recent analysis)
# We use 'priority=['aws']' to ensure it pulls from the Amazon cloud
try:
    H = Herbie(model='hrrr', product='sfc', fxx=0, priority=['aws'])
    ds = H.xarray("TMP:2 m") # Pulls 2-meter temperature
    
    # 2. Process Data (Kelvin to Fahrenheit)
    temp_f = (ds.t2m - 273.15) * 9/5 + 32

    # 3. Create the Visualization
    fig = plt.figure(figsize=(12, 7))
    ax = plt.axes(projection=ds.herbie.crs)
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    
    # Add data to map
    plot = ax.pcolormesh(ds.longitude, ds.latitude, temp_f, 
                         transform=ccrs.PlateCarree(), cmap='RdYlBu_r')
    
    plt.colorbar(plot, label="Temperature (°F)", orientation='vertical', pad=0.02)
    plt.title(f"HRRR 2m Temperature - Valid: {H.date}")

    # 4. Save the image
    plt.savefig("latest_hrrr.png", dpi=150, bbox_inches='tight')
    print("Map generated successfully.")

except Exception as e:
    print(f"Error: {e}")
