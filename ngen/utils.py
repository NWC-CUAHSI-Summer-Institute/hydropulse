import geopandas as gpd
import fiona
import folium
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from IPython.display import Markdown, display
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os

def clip_hydrofabric(hydrofabric_gpkg, boundary_gpkg, output_dir, name):
    """
    Clip all layers in the hydrofabric GPKG using the boundary from another GPKG layer,
    and save the result to a new GPKG in the output directory.

    Parameters:
        hydrofabric_gpkg (str): Path to the input hydrofabric GPKG file.
        boundary_gpkg (str): Path to the boundary GPKG file (must have only one layer).
        output_dir (str): Output directory where the result GPKG will be saved.
        name (str): Output file name (without .gpkg).
    """
    os.makedirs(output_dir, exist_ok=True)
    output_gpkg = os.path.join(output_dir, f"{name}.gpkg")

    # Automatically detect the only boundary layer
    boundary_layer = fiona.listlayers(boundary_gpkg)[0]
    boundary = gpd.read_file(boundary_gpkg, layer=boundary_layer)
    
    # Reproject boundary to EPSG:5070 (NAD83 / Conus Albers)
    boundary = boundary.to_crs(epsg=5070)

    # List all layers in the hydrofabric GPKG
    layers = fiona.listlayers(hydrofabric_gpkg)

    for layer in layers:
        print(f"Clipping layer: {layer}")
        gdf = gpd.read_file(hydrofabric_gpkg, layer=layer)
    
        # Skip non-spatial or corrupted layers
        if not hasattr(gdf, 'geometry') or gdf.geometry.name not in gdf.columns:
            print(f"Skipping non-spatial or corrupted layer: {layer}")
            continue
    
        gdf = gdf.to_crs(boundary.crs)
    
        try:
            clipped = gpd.overlay(gdf, boundary, how='intersection')
        except Exception as e:
            print(f"Error clipping {layer}: {e}")
            continue
    
        if not clipped.empty:
            clipped.to_file(output_gpkg, layer=layer, driver="GPKG")
        else:
            print(f"Layer {layer} has no overlapping features after clipping.")

#Plot the ngen hydrofabric
def plot_ngen_hydrofabric(hf_filepath, catchment_id=None, gage_id=None, output_dir="../outputs"):
    """
    Plot the NGEN HydroFabric elements for a given catchment or gage.

    Parameters:
    - divides: GeoDataFrame of catchment divides
    - flowpaths: GeoDataFrame of flowpaths (reaches)
    - nexus: GeoDataFrame of nexus points
    - catchment_id: str, optional catchment identifier (e.g., "2430598")
    - gage_id: str, optional gage identifier (e.g., "08074000")
    - output_dir: str, path to directory to save output image
    """
    if not (catchment_id or gage_id):
        raise ValueError("Must provide either a catchment_id or gage_id.")
    if catchment_id and gage_id:
        raise ValueError("Provide only one of catchment_id or gage_id, not both.")

    id_str = catchment_id if catchment_id else gage_id
    
    title_type = "catchment" if catchment_id else "gage"
    
    divides = gpd.read_file(hf_filepath, layer="divides")
    flowpaths = gpd.read_file(hf_filepath, layer="flowpaths")
    nexus = gpd.read_file(hf_filepath, layer="nexus")

    fig, ax = plt.subplots(figsize=(8, 8))

    divides.plot(ax=ax, color='green', alpha=0.3, edgecolor='darkgreen')
    flowpaths.plot(ax=ax, color='blue', edgecolor='blue', alpha=0.5, label='HydroFabric Reaches')
    nexus.plot(ax=ax, color='red', edgecolor='red', alpha=0.5, label='HydroFabric Nexus')

    ax.set_title(f'NGEN HydroFabric at {title_type} {id_str}')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    # Rotate y-axis tick labels
    for label in ax.get_yticklabels():
        label.set_rotation(90)

    ax.legend()
    output_path = f"{output_dir}/Ngen_HF_{id_str}.png"
    plt.savefig(output_path, dpi=500, bbox_inches='tight')
    plt.show()

#Plot the hydrifabric
def plot_hydrofabric(hf_filepath):
    """
    Load and display an interactive Folium map from a HydroFabric GPKG file.
    Layers used: divides, flowpaths, nexus

    Parameters:
        hf_filepath (str): Path to the HydroFabric GPKG file.

    Returns:
        folium.Map object
    """

    def clean_gdf_for_folium(gdf):
        gdf = gdf.copy()
        for col in gdf.columns:
            if gdf[col].dtype.name.startswith("datetime64") or isinstance(gdf[col].iloc[0], pd.Timestamp):
                gdf[col] = gdf[col].astype(str)
        return gdf

    def valid_fields(gdf):
        return [col for col in gdf.columns if col != 'geometry']

    # Load & reproject
    divides = gpd.read_file(hf_filepath, layer="divides").to_crs(epsg=4326)
    flowpaths = gpd.read_file(hf_filepath, layer="flowpaths").to_crs(epsg=4326)
    nexus = gpd.read_file(hf_filepath, layer="nexus").to_crs(epsg=4326)

    # Clean datetime fields
    divides = clean_gdf_for_folium(divides)
    flowpaths = clean_gdf_for_folium(flowpaths)
    nexus = clean_gdf_for_folium(nexus)

    # Center map on bounding box center of divides
    bounds = divides.total_bounds  # (minx, miny, maxx, maxy)
    center = [(bounds[1] + bounds[3]) / 2, (bounds[0] + bounds[2]) / 2]

    # Create base map
    m = folium.Map(location=center, zoom_start=11, tiles="CartoDB positron")

    # Add Divides
    folium.GeoJson(
        divides,
        name='Divides',
        style_function=lambda x: {
            'fillColor': 'green', 'color': 'darkgreen', 'weight': 1, 'fillOpacity': 0.3
        },
        highlight_function=lambda x: {
            'fillColor': '#ffffcc', 'color': 'orange', 'weight': 3, 'fillOpacity': 0.7
        },
        popup=folium.GeoJsonPopup(
            fields=valid_fields(divides),
            aliases=valid_fields(divides),
            labels=True,
            sticky=True,
            max_width=300
        )
    ).add_to(m)

    # Add Flowpaths
    folium.GeoJson(
        flowpaths,
        name='Flowpaths',
        style_function=lambda x: {
            'color': 'blue', 'weight': 2, 'opacity': 0.5
        },
        highlight_function=lambda x: {
            'color': 'red', 'weight': 3, 'opacity': 0.9
        },
        popup=folium.GeoJsonPopup(
            fields=valid_fields(flowpaths),
            aliases=valid_fields(flowpaths),
            labels=True,
            sticky=True,
            max_width=300
        )
    ).add_to(m)

    # Add Nexus points
    for _, row in nexus.iterrows():
        folium.CircleMarker(
            location=[row.geometry.y, row.geometry.x],
            radius=5,
            color='red',
            fill=True,
            fill_color='red',
            fill_opacity=0.8,
            popup=folium.Popup(
                html="<br>".join([
                    f"<b>{col}:</b> {row[col]}" for col in valid_fields(nexus)
                ]),
                max_width=300
            )
        ).add_to(m)

    folium.LayerControl().add_to(m)
    return m

#Get the subset of the hydrofabric
def run_ngiab_subset(gage_id=None, catchment_id=None):
    """
    Run ngiab_data_cli subset command for gage or catchment.

    Parameters:
        gage_id (str): USGS gage site number (e.g., '02464000')
        catchment_id (str): Catchment ID (e.g., '2430598')
    """
    if gage_id and catchment_id:
        raise ValueError("Please provide either gage_id or catchment_id, not both.")

    if gage_id:
        cmd = ["python", "-m", "ngiab_data_cli", "-i", f"gage-{gage_id}", "--subset"]
    elif catchment_id:
        cmd = ["python", "-m", "ngiab_data_cli", "-i", f"cat-{catchment_id}", "--subset"]
    else:
        raise ValueError("You must specify either a gage_id or a catchment_id.")

    print("Running command:", " ".join(cmd))
    subprocess.run(cmd, check=True)

#Get forcing
def get_ngiab_forcings(gage_id=None, catchment_id=None, cat_ids = None, start=None, end=None):
    """
    Run ngiab_data_cli for either gage or catchment with optional date range and full flag.

    Parameters:
        gage_id (str): USGS gage site number (e.g., '02464000')
        catchment_id (str): Catchment ID (e.g., '2430598')
        start (str): Start date in YYYY-MM-DD format
        end (str): End date in YYYY-MM-DD format
        full (bool): If True, include the -f (full) flag
        cat_ids : Is for multiple catchments (only particular one)
    """
    if sum([gage_id is not None, catchment_id is not None, bool(cat_ids)]) > 1:
        raise ValueError("Provide only one of: gage_id, catchment_id, or cat_ids.")
    if not any([gage_id is not None, catchment_id is not None, bool(cat_ids)]):
        raise ValueError("You must specify exactly one of: gage_id, catchment_id, or cat_ids.")

    # Base ID
    id_str = f"gage-{gage_id}" if gage_id else (f"cat-{catchment_id}" if catchment_id else f"wb-{cat_ids[0]}")

    # Start building command
    cmd = ["python", "-m", "ngiab_data_cli", "-i", id_str]

    cmd.append("-f")
    # Add date range if provided
    if start:
        cmd += ["--start", start]
    if end:
        cmd += ["--end", end]

    print("Running command:", " ".join(cmd))
    subprocess.run(cmd, check=True)

#Preparing the configuration files, 
def run_ngiab_realization(gage_id=None, catchment_id=None, cat_ids = None, start=None, end=None):
    """
    Run ngiab_data_cli with -r and show the config directory structure cleanly.

    Parameters:
        gage_id (str): USGS gage site number
        catchment_id (str): Catchment ID
        start (str): Start date (YYYY-MM-DD)
        end (str): End date (YYYY-MM-DD)
    """
    if sum([gage_id is not None, catchment_id is not None, bool(cat_ids)]) > 1:
        raise ValueError("Provide only one of: gage_id, catchment_id, or cat_ids.")
    if not any([gage_id is not None, catchment_id is not None, bool(cat_ids)]):
        raise ValueError("You must specify exactly one of: gage_id, catchment_id, or cat_ids.")

    id_str = f"gage-{gage_id}" if gage_id else (f"cat-{catchment_id}" if catchment_id else f"wb-{cat_ids[0]}")
    
    output_base = "/home/jovyan/ngiab_preprocess_output"
    config_dir = os.path.join(output_base, id_str, "config")

    cmd = ["python", "-m", "ngiab_data_cli", "-i", id_str, "-r"]
    if start:
        cmd += ["--start", start]
    if end:
        cmd += ["--end", end]

    print("Running command:", " ".join(cmd))
    subprocess.run(cmd, check=True)

#Plotting the CFE results 
def plot_CFE_results(
    outletID_value,
    cat_id,
    id_type, 
    start_time,
    end_time,
    variables,
    output_dir,
    base_dir='../ngiab_preprocess_output'
):
    """
    Plot selected NGEN variables over time and save the plot using the catchment or gage ID.
    """
    id_str = f"{id_type}-{outletID_value}"
    ID = f"cat-{cat_id}"
    csv_path = Path(base_dir) / id_str / 'outputs' / 'ngen' / f"{ID}.csv"

    if not csv_path.exists():
        print(f"CSV not found: {csv_path}")
        return

    # Load and parse time
    df = pd.read_csv(csv_path, parse_dates=['Time'])
    df['Time'] = pd.to_datetime(df['Time'])

    # Filter by time
    start_time = pd.to_datetime(start_time)
    end_time = pd.to_datetime(end_time)
    year_label = start_time.year

    df_filtered = df[(df['Time'] >= start_time) & (df['Time'] <= end_time)]

    if df_filtered.empty:
        print(f"No data in selected time range for {id_str}")
        return

    # Plotting
    plt.figure(figsize=(7, 3))
    handles = []
    labels = []
    
    for var in variables:
        if var in df_filtered.columns:
            if var == 'GIUH_RUNOFF':
                line, = plt.plot(
                    df_filtered['Time'], df_filtered[var],
                    linestyle='-.', linewidth=2.2
                )
            else:
                line, = plt.plot(
                    df_filtered['Time'], df_filtered[var],
                    linestyle='-', linewidth=1.2
                )
            handles.append(line)
            labels.append(var)
        else:
            print(f"Variable '{var}' not found in {ID}.csv")
    
    # Axis formatting
    ax = plt.gca()
    tick_locs = np.linspace(df_filtered['Time'].min().value, df_filtered['Time'].max().value, 5)
    ax.set_xticks(pd.to_datetime(tick_locs))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d-%m'))
    
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel(f'Date (Year: {year_label})', fontsize=12)
    plt.ylabel('Depth (in m)', fontsize=12)
    
    # Add cat ID as first legend entry
    handles.insert(0, plt.Line2D([], [], color='none'))
    labels.insert(0, f"NextGen Catchment {cat_id}")
    
    plt.legend(handles, labels, fontsize=10, frameon=False)
    
    plt.grid(True, linestyle='-.', alpha=0.7)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    
    plt.tight_layout()
    
    # Save plot
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    save_path = output_dir / f"{ID}_timeseries.png"
    plt.savefig(save_path, dpi=500)
    print(f"Plot saved: {save_path}")
    plt.show()
