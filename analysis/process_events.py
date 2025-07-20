import pandas as pd
from pathlib import Path
from calculate_area import calculate_area
from select_giuh_timeline import select_giuh_timeline
from calculate_indices import calculate_indices
from scipy.stats import genextreme as gev
import numpy as np

def extract_annual_max_rainfall(catchment_id, cfe_results_dir):
    """Extract annual maximum 1-hour rainfall from CFE CSV for a catchment from 2005 to 2021."""
    try:
        csv_file = Path(cfe_results_dir) / f"{catchment_id}.csv"
        if not csv_file.exists():
            print(f"[SKIPPED] CFE file not found for catchment {catchment_id}")
            return None

        # Load CFE CSV
        df = pd.read_csv(csv_file, parse_dates=["Time"])
        df["Time"] = pd.to_datetime(df["Time"], errors="coerce")
        df = df.dropna(subset=["Time"])
        df["RAIN_RATE"] = pd.to_numeric(df["RAIN_RATE"], errors="coerce")
        
        # Filter data from 2005 to 2021
        df = df[(df["Time"].dt.year >= 2005) & (df["Time"].dt.year <= 2021)]
        if df.empty:
            print(f"[SKIPPED] No data from 2005 to 2021 for catchment {catchment_id}")
            return None

        # Group by year and find maximum RAIN_RATE
        annual_max_rainfall = df.groupby(df["Time"].dt.year)["RAIN_RATE"].max().dropna()
        if len(annual_max_rainfall) < 10:
            print(f"[SKIPPED] Fewer than 10 years of data for catchment {catchment_id}")
            return None

        return annual_max_rainfall

    except Exception as e:
        print(f"[ERROR] Failed to extract annual max rainfall for catchment {catchment_id}: {e}")
        return None

def calculate_return_period(catchment_id, rainfall_value, cfe_results_dir):
    """Calculate return period for a given rainfall value using GEV distribution."""
    try:
        # Extract annual maximum rainfall
        annual_max_rainfall = extract_annual_max_rainfall(catchment_id, cfe_results_dir)
        if annual_max_rainfall is None:
            return np.nan

        # Fit GEV distribution
        params = gev.fit(annual_max_rainfall)
        c, loc, scale = params  # Shape, location, scale

        # Calculate exceedance probability
        exceedance_prob = 1 - gev.cdf(rainfall_value, c, loc, scale)
        
        # Avoid division by zero
        if exceedance_prob <= 0:
            print(f"[WARNING] Exceedance probability is zero for catchment {catchment_id}, rainfall {rainfall_value}")
            return np.nan

        # Calculate return period
        return_period = 1 / exceedance_prob
        return return_period

    except Exception as e:
        print(f"[ERROR] Failed to calculate return period for catchment {catchment_id}: {e}")
        return np.nan

def process_events(events_file, catchments_dir, cfe_results_dir, output_file):
    """Process all events for unique VALID2 dates and generate output CSV, including return period."""
    # Read events CSV
    events_df = pd.read_csv(events_file, parse_dates=["VALID2"])
    
    # Validate required columns
    required_columns = {"nextgen_catchment", "VALID2", "damagePercent"}
    if not required_columns.issubset(events_df.columns):
        missing = required_columns - set(events_df.columns)
        raise ValueError(f"Missing required columns in Events.csv: {missing}")
    
    results = []
    log_file_path = Path("skipped_catchments.log")
    
    # Group by unique VALID2 dates
    for valid_date in events_df["VALID2"].unique():
        # Round VALID2 to the nearest next hour
        valid_date = pd.to_datetime(valid_date)
        updated_date = valid_date.ceil("h")
        
        # Get all rows for this VALID2 date
        date_rows = events_df[events_df["VALID2"] == valid_date]
        
        for _, row in date_rows.iterrows():
            catchment_id = row["nextgen_catchment"]
            damage_percent = row["damagePercent"]
            
            try:
                # Calculate area
                area_sqkm = calculate_area(catchment_id, catchments_dir)
                
                # Select GIUH timeline and get alpha1, rainrate, giuh_runoff, and starttime
                result = select_giuh_timeline(catchment_id, updated_date, cfe_results_dir)
                if result is None:
                    error_msg = f"Skipping catchment {catchment_id} for date {valid_date} due to error in select_giuh_timeline"
                    print(error_msg)
                    with open(log_file_path, "a") as log_file:
                        log_file.write(f"{error_msg}\n")
                    continue
                window_df, alpha1, rainrate, giuh_runoff, starttime = result
                
                # Calculate indices
                runoff_volume, rainfall_volume, pfi_1, pfi_2 = calculate_indices(window_df, area_sqkm)
                
                # Calculate time difference between updated_event_date and starttime (in hours)
                time_diff_hours = (updated_date - starttime).total_seconds() / 3600
                
                # Calculate return period for rainrate
                return_period = calculate_return_period(catchment_id, rainrate, cfe_results_dir)
                
                # Store results
                result = {
                    "Event_date": valid_date,
                    "updated_event_date": updated_date,
                    "catchment": catchment_id,
                    "area_sqkm": area_sqkm,
                    "runoff_volume": runoff_volume,
                    "rainfall_volume": rainfall_volume,
                    "PFI_1": pfi_1,
                    "PFI_2": pfi_2,
                    "damage": damage_percent,
                    "rainrate": rainrate,
                    "giuh_runoff": giuh_runoff,
                    "alpha1": alpha1,
                    "starttime": starttime,
                    "time_diff_hours": time_diff_hours,
                    "rainfall_return_period_years": return_period
                }
                results.append(result)
                
            except (FileNotFoundError, ValueError) as e:
                error_msg = f"Skipping catchment {catchment_id} for date {valid_date} due to error: {e}"
                print(error_msg)
                with open(log_file_path, "a") as log_file:
                    log_file.write(f"{error_msg}\n")
                continue
    
    # Create results DataFrame and save to CSV
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
    return results_df