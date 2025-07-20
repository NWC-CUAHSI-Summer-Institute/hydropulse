import pandas as pd
import numpy as np

def calculate_indices(window_df, area_sqkm):
    """Calculate PFI-1 (runoff_volume / rainfall_volume) and PFI-2 (peak GIUH / peak Rainfall)."""
    # Check if window_df is empty or contains only NaN values
    if window_df.empty or window_df["GIUH_RUNOFF"].isna().all() or window_df["RAIN_RATE"].isna().all():
        raise ValueError("Invalid window_df: empty or contains only NaN values for GIUH_RUNOFF or RAIN_RATE")
    
    # Count hours where GIUH_RUNOFF is not zero (optional for analysis)
    non_zero_hours = len(window_df[window_df["GIUH_RUNOFF"] > 0])
    
    # Calculate total runoff and rainfall volumes
    runoff_volume = window_df["GIUH_RUNOFF"].sum()
    rainfall_volume = window_df["RAIN_RATE"].sum()

    # Compute PFI-1, avoiding division by zero
    pfi_1 = runoff_volume / rainfall_volume if rainfall_volume > 0 else 0.0

    # Compute PFI-2, avoiding division by zero
    peak_giuh = window_df["GIUH_RUNOFF"].max()
    peak_rainfall = window_df["RAIN_RATE"].max()
    pfi_2 = peak_giuh / peak_rainfall if peak_rainfall > 0 else 0.0

    return runoff_volume, rainfall_volume, pfi_1, pfi_2
