from skyfield.api import Topos, load
from skyfield.sgp4lib import EarthSatellite

# TLE data
line1 = '1 39215U 13038A   23283.17449397  .00000133  00000+0  00000+0 0  9992'
line2 = '2 39215   2.8519   9.5912 0003790 204.4114 252.2040  1.00272111 37435'

# Create satellite object from the TLE lines
satellite = EarthSatellite(line1, line2)

# Load a timescale and ask for the current time
ts = load.timescale()
t = ts.now()

# Get the geocentric position of the satellite
geocentric = satellite.at(t)

# Extract sub-satellite point (latitude and longitude)
subpoint = geocentric.subpoint()
lat_deg = subpoint.latitude.degrees
lon_deg = subpoint.longitude.degrees

print(f"Latitude: {lat_deg}°")
print(f"Longitude: {lon_deg}°")
