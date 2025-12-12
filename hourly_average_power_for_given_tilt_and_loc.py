#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# --- Constants and Location ---
lat_deg = 46.5                # Latitude for Lausanne
lat = np.deg2rad(lat_deg)
tilt_deg = 20                 # Fixed panel tilt (degrees from horizontal)
tilt = np.deg2rad(tilt_deg)
ridge_orientation = 0         # North-south oriented ridge (so panels face east and west)

# For a symmetric two-pan roof:
#   Panel 1 normal: ridge_orientation - 90°  -> east-facing
#   Panel 2 normal: ridge_orientation - 270° -> west-facing (equivalent to +90°)
panel1_az_deg = ridge_orientation - 90  
panel2_az_deg = ridge_orientation - 270 
panel1_az = np.deg2rad(panel1_az_deg)
panel2_az = np.deg2rad(panel2_az_deg)

# --- Prepare to record hourly instantaneous power ---
times = []       # Time in fractional day-of-year (e.g., 5.5 means day 5 at 12:00)
power_inst = []  # Instantaneous power output (W/m²) for each hour

# --- Loop over each day and hour of the year ---
for day in range(1, 366):
    # Compute solar declination (δ) in radians
    delta = np.deg2rad(23.45) * np.sin(2 * np.pi * (284 + day) / 365)
    for hour in range(24):
        # Hour angle ω (in radians): ω = (hour - 12)*15°
        omega = np.deg2rad((hour - 12) * 15)
        
        # Solar altitude (α): sin(α) = sin(lat)*sin(δ) + cos(lat)*cos(δ)*cos(ω)
        sin_alpha = np.sin(lat) * np.sin(delta) + np.cos(lat) * np.cos(delta) * np.cos(omega)
        if sin_alpha <= 0:
            # Sun is below horizon, so no power generated
            power = 0
        else:
            alpha = np.arcsin(sin_alpha)
            # Solar zenith angle (θ_z)
            theta_z = np.pi / 2 - alpha

            # Solar azimuth calculation:
            cos_az = (np.sin(delta) - np.sin(alpha) * np.sin(lat)) / (np.cos(alpha) * np.cos(lat))
            cos_az = np.clip(cos_az, -1, 1)
            az = np.arccos(cos_az)
            if omega > 0:
                az = 2 * np.pi - az
            # Convert solar azimuth to “relative to south” by subtracting π.
            solar_az_rel = az - np.pi

            # --- Direct Beam Irradiance Model ---
            am = 1 / sin_alpha  # Air mass approximation
            I_beam = 1000 * (0.7 ** (am ** 0.678))  # Clear-sky beam irradiance (W/m²)
            
            # --- Incidence Angle on Each Panel ---
            cos_inc1 = np.cos(theta_z) * np.cos(tilt) + np.sin(theta_z) * np.sin(tilt) * np.cos(solar_az_rel - panel1_az)
            cos_inc2 = np.cos(theta_z) * np.cos(tilt) + np.sin(theta_z) * np.sin(tilt) * np.cos(solar_az_rel - panel2_az)
            cos_inc1 = max(cos_inc1, 0)
            cos_inc2 = max(cos_inc2, 0)
            
            # Instantaneous power (W/m²) from both panels
            power = I_beam * (cos_inc1 + cos_inc2)
        
        # Record the time (day + hour/24) and corresponding instantaneous power
        times.append(day + hour / 24.0)
        power_inst.append(power)

# Convert lists to numpy arrays for plotting
times = np.array(times)
power_inst = np.array(power_inst)

# --- Plotting the Hourly Instantaneous Power ---
plt.figure(figsize=(12, 6))
plt.plot(times, power_inst, lw=0.5, label="Instantaneous Power")
plt.xlabel("Day of Year")
plt.ylabel("Instantaneous Power (W/m²)")
plt.title("Hourly Instantaneous Power Output per m² in Lausanne\nTilt: {}°, Ridge Orientation: North-South".format(tilt_deg))
plt.grid(True)
plt.legend()
plt.savefig("hourly_power_plot.png")
plt.show()
