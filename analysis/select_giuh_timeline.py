import pandas as pd
from pathlib import Path
from datetime import timedelta
import numpy as np

def select_giuh_timeline(catchment_id, valid_date, cfe_results_dir):
    """Select GIUH_RUNOFF for 18h from starttime where alpha1 >= 11.2. Return window_df, alpha1, rainrate, giuh_runoff, and starttime."""
    csv_file = Path(cfe_results_dir) / f"{catchment_id}.csv"
    if not csv_file.exists():
        print(f"[SKIPPED] File not found for catchment {catchment_id}")
        return None

    try:
        df = pd.read_csv(csv_file, parse_dates=["Time"])
        df["Time"] = pd.to_datetime(df["Time"], errors="coerce")
        df = df.set_index("Time").sort_index()
        df.index = df.index.tz_localize(None)
        df["GIUH_RUNOFF"] = pd.to_numeric(df["GIUH_RUNOFF"], errors="coerce")
        df["RAIN_RATE"] = pd.to_numeric(df["RAIN_RATE"], errors="coerce")
    except Exception as e:
        print(f"[ERROR] Failed to load or process file for catchment {catchment_id}: {e}")
        return None

    if not pd.api.types.is_datetime64_any_dtype(valid_date):
        valid_date = pd.to_datetime(valid_date, errors="coerce")
        if pd.isna(valid_date):
            print(f"[SKIPPED] Invalid date for catchment {catchment_id}")
            return None

    start_time = valid_date - timedelta(hours=24)

    try:
        search_df = df.loc[start_time:valid_date]
    except Exception as e:
        print(f"[ERROR] KeyError while slicing data for catchment {catchment_id} from {start_time} to {valid_date}: {e}")
        return None

    if search_df.empty:
        print(f"[SKIPPED] No GIUH_RUNOFF data for catchment {catchment_id} from {start_time} to {valid_date}")
        return None

    starttime = None
    alpha1 = np.nan
    rainrate = np.nan
    giuh_runoff = np.nan

    for i in range(len(search_df)):
        current_giuh = search_df["GIUH_RUNOFF"].iloc[i]
        current_rainrate = search_df["RAIN_RATE"].iloc[i]
        if not pd.isna(current_giuh) and not pd.isna(current_rainrate) and current_rainrate > 0:
            current_alpha1 = (current_giuh / current_rainrate) * 100
            if current_alpha1 >= 11.2:
                starttime = search_df.index[i]
                giuh_runoff = current_giuh
                rainrate = current_rainrate
                alpha1 = current_alpha1
                break

    if starttime is None:
        print(f"[SKIPPED] No alpha1 >= 11.2 found for catchment {catchment_id} from {start_time} to {valid_date}")
        return None

    end_time = starttime + timedelta(hours=18)

    try:
        window_df = df.loc[starttime:end_time]
    except Exception as e:
        print(f"[ERROR] Failed to get window from {starttime} to {end_time} for {catchment_id}: {e}")
        return None

    if window_df.empty:
        print(f"[SKIPPED] No GIUH_RUNOFF in 18-hour window from {starttime} to {end_time} for {catchment_id}")
        return None

    return window_df, alpha1, rainrate, giuh_runoff, starttime