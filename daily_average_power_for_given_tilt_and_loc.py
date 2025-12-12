#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# --- Constants and Location ---
lat_deg = 46.5                # Latitude for Lausanne
lat = np.deg2rad(lat_deg)
tilt_deg = 20                 # Fixed panel tilt (degrees from horizontal)
tilt = np.deg2rad(tilt_deg)
ridge_orientation = 0         # North-south oriented ridge (0° from north)

# For a symmetric two-pan roof:
#   Panel 1 normal: ridge_orientation - 90° (relative-to-south)
#   Panel 2 normal: ridge_orientation - 270° (or equivalently, +90°)
panel1_az_deg = ridge_orientation - 90  # = -90° -> effectively facing east
panel2_az_deg = ridge_orientation - 270 # = -270° which is equivalent to 90° -> facing west
panel1_az = np.deg2rad(panel1_az_deg)
panel2_az = np.deg2rad(panel2_az_deg)

# --- Time Setup ---
days = np.arange(1, 366)             # Days 1 to 365
daily_average_power = np.zeros(len(days))  # To store the daily average (W/m²)

# --- Loop over each day ---
for idx, day in enumerate(days):
    daily_energy = 0.0  # Wh/m² collected during the day
    for hour in range(24):
        # --- Solar Position Calculations ---
        # Solar declination δ (in radians)
        delta = np.deg2rad(23.45) * np.sin(2 * np.pi * (284 + day) / 365)
        # Hour angle ω (in radians): ω = (hour - 12)*15°
        omega = np.deg2rad((hour - 12) * 15)
        
        # Solar altitude α: sin(α) = sin(lat)*sin(δ) + cos(lat)*cos(δ)*cos(ω)
        sin_alpha = np.sin(lat) * np.sin(delta) + np.cos(lat) * np.cos(delta) * np.cos(omega)
        # Skip hours when sun is below horizon
        if sin_alpha <= 0:
            continue
        alpha = np.arcsin(sin_alpha)
        # Solar zenith angle θ_z
        theta_z = np.pi/2 - alpha

        # Solar azimuth calculation:
        # cos(azimuth) = (sin(δ) - sin(α)*sin(lat)) / (cos(α)*cos(lat))
        cos_az = (np.sin(delta) - np.sin(alpha)*np.sin(lat)) / (np.cos(alpha)*np.cos(lat))
        cos_az = np.clip(cos_az, -1, 1)
        az = np.arccos(cos_az)
        # Adjust for hour angle (sun in the west when ω > 0)
        if omega > 0:
            az = 2 * np.pi - az
        # Convert solar azimuth to “relative to south” by subtracting π.
        solar_az_rel = az - np.pi

        # --- Direct Beam Irradiance Model ---
        # Air mass approximation: AM = 1/sin(α)
        am = 1 / sin_alpha
        # Clear-sky beam irradiance (W/m²)
        I_beam = 1000 * (0.7 ** (am ** 0.678))
        
        # --- Incidence Angle on Each Panel ---
        # For panel 1:
        cos_inc1 = np.cos(theta_z)*np.cos(tilt) + np.sin(theta_z)*np.sin(tilt)*np.cos(solar_az_rel - panel1_az)
        cos_inc1 = max(cos_inc1, 0)
        # For panel 2:
        cos_inc2 = np.cos(theta_z)*np.cos(tilt) + np.sin(theta_z)*np.sin(tilt)*np.cos(solar_az_rel - panel2_az)
        cos_inc2 = max(cos_inc2, 0)
        
        # Total instantaneous power (W/m²) is the sum from both panels
        I_total = I_beam * (cos_inc1 + cos_inc2)
        # Add the energy contribution for this hour (Wh/m²)
        daily_energy += I_total

    # Compute daily average hourly power (divide daily energy by 24 hours)
    daily_average_power[idx] = daily_energy / 24

# --- Plotting ---
plt.figure(figsize=(10, 6))
plt.plot(days, daily_average_power, label=f"Tilt = {tilt_deg}°")
plt.xlabel("Day of Year")
plt.ylabel("Average Hourly Power (W/m²)")
plt.title("Daily Average Hourly Power Output per m² (Clear-Sky Model) in Lausanne")
plt.grid(True)
plt.legend()
plt.savefig("daily_power_plot.png")
plt.show()
