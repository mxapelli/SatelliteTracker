from skyfield.api import Topos, load
from skyfield.sgp4lib import EarthSatellite
import datetime

def compute_satellite_positions(line1, line2):
    # Create satellite object from the TLE lines
    satellite = EarthSatellite(line1, line2)

    # Load a timescale
    ts = load.timescale()
    
    # Start from the current time
    t0 = ts.now()

    # Calculate the end time (1 hour from now)
    end_time = t0.utc_datetime() + datetime.timedelta(hours=1)
    t_end = ts.utc(end_time.year, end_time.month, end_time.day, end_time.hour, end_time.minute, end_time.second)

    positions = []

    # Loop every 5 seconds from now until 1 hour
    interval = datetime.timedelta(seconds=5)
    while t0.utc_datetime() <= t_end.utc_datetime():
        geocentric = satellite.at(t0)
        subpoint = geocentric.subpoint()
        lat_deg = subpoint.latitude.degrees
        lon_deg = subpoint.longitude.degrees
        
        positions.append((t0.utc_datetime(), lat_deg, lon_deg))

        # Increment the time by 5 seconds
        t0 = ts.utc((t0.utc_datetime() + interval).year, 
                    (t0.utc_datetime() + interval).month, 
                    (t0.utc_datetime() + interval).day, 
                    (t0.utc_datetime() + interval).hour, 
                    (t0.utc_datetime() + interval).minute, 
                    (t0.utc_datetime() + interval).second)
    
    return positions


