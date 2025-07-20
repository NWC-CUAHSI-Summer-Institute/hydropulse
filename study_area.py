import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from pathlib import Path
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import warnings
warnings.filterwarnings("ignore")

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def plot_LSRincidents(
    gpkg_path,  
    csv_path, 
    layer_name='divides',
    divide_id_col='divide_id', 
    catchment_col='nextgen_catchment',
    damage_col='damagePercent',
    crs_proj='EPSG:5070'
):
    # Load catchment polygons
    gdf = gpd.read_file(gpkg_path, layer=layer_name)

    if divide_id_col not in gdf.columns:
        raise ValueError(f"'{divide_id_col}' not found in GPKG.")

    # Project catchments to projected CRS
    gdf_proj = gdf.to_crs(crs_proj)
    gdf_proj['centroid'] = gdf_proj.geometry.centroid
    centroid_gdf = gdf_proj.set_geometry('centroid')

    # Load CSV and compute unique damage counts per catchment
    df = pd.read_csv(csv_path).dropna(subset=[damage_col])
    unique_counts = (
        df.groupby(catchment_col)[damage_col]
          .nunique()
          .reset_index()
          .rename(columns={catchment_col: divide_id_col, damage_col: "unique_damage_count"})
    )

    merged = centroid_gdf.merge(unique_counts, on=divide_id_col, how='left')
    merged['unique_damage_count'] = merged['unique_damage_count'].fillna(0)

    # Setup color scale
    vmin = merged['unique_damage_count'].min()
    vmax = merged['unique_damage_count'].max()
    print(vmin, vmax)
    bounds = np.linspace(vmin, vmax, 4)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # Convert to lat/lon for plotting with Cartopy
    merged_latlon = merged.to_crs("EPSG:4326")

    # Plot with Cartopy
    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=ccrs.LambertConformal())
    ax.set_extent([-120, -74, 24, 50], crs=ccrs.PlateCarree())

    # Add base features
    ax.add_feature(cfeature.BORDERS, linewidth=1)
    ax.add_feature(cfeature.STATES, linewidth=0.2)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)

    # Plot centroids
    sc = ax.scatter(
        merged_latlon.geometry.x,
        merged_latlon.geometry.y,
        c=merged_latlon['unique_damage_count'],
        cmap='coolwarm',
        norm=norm,
        s=100,
        alpha=0.8,
        edgecolors='black',
        linewidths=0.3,
        transform=ccrs.PlateCarree()
    )

    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, orientation='vertical', shrink=0.4, pad=0.01, aspect=20)
    cbar.set_label('Number of LSR events with NFIP damage \nrecords per NextGen catchment', fontsize=12)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_ticks(bounds)
    cbar.set_ticklabels([f"{int(b)}" for b in bounds])
    area_data = merged_latlon['areasqkm']
    bins = [0, 5, 10, 15, np.inf]
    labels = ['<5', '5–10', '10–15', '>15']
    counts, _ = np.histogram(area_data, bins=bins)

    inset_ax = inset_axes(
        ax,
        width="18%", height="25%",
        loc='lower left',
        bbox_to_anchor=(0.08, 0.08, 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=0
    )

    inset_ax.patch.set_alpha(0.5)

    bar_width = 0.1
    positions = np.arange(len(labels)) * bar_width  

    bar_container = inset_ax.bar(
        positions, counts, width=bar_width, align='edge',
        color='azure', edgecolor='black', linewidth=1.2
    )

    # Show only left and bottom spines
    for spine_pos, spine in inset_ax.spines.items():
        spine.set_visible(spine_pos in ['left', 'bottom'])

    # Tight label placement
    inset_ax.set_xlabel('Area Range (km²)', fontsize=11, labelpad=1)
    inset_ax.xaxis.label.set_bbox(dict(facecolor='white', edgecolor='none', alpha=0.8, boxstyle='round,pad=0'))
    
    inset_ax.set_ylabel('Count', fontsize=11, labelpad=0.5)
    inset_ax.yaxis.label.set_bbox(dict(facecolor='white', edgecolor='none', alpha=0.8, boxstyle='round,pad=0'))


    # Y-axis ticks: bottom, mid, top
    y_max = max(counts)
    inset_ax.set_yticks([0, y_max // 2, y_max])
    inset_ax.set_yticklabels([0, y_max // 2, y_max], fontsize=11)

    # No x-ticks (we label inside bars)
    inset_ax.set_xticks([])

    # Area range labels inside each bar, vertically
    for i, bar in enumerate(bar_container):
        height = bar.get_height()
        if height > 0:
            inset_ax.text(
                bar.get_x() + bar.get_width() / 2, height * 0.5, labels[i],
                ha='center', va='center', fontsize=10, rotation=90, color='black'
            )



    # Summary text
    total_points = merged_latlon.shape[0]
    info_text = f"NextGen Catchments: {total_points}"
    ax.text(
        0.98, 0.98, info_text,
        transform=ax.transAxes,
        ha='right', va='top',
        fontsize=11,
        bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, boxstyle='round,pad=0.3')
    )

    plt.tight_layout()
    plt.savefig("../outputs/study_area.png", dpi=500)
    plt.show()

#These are the selected Nextgen catchments for the study area on basis of the number of LSR events
ngen_cat = Path("../data/FinalData/Nextgen_cat_subset.gpkg")

#These are the LSR events with NFIP damage records
lsr_events = Path("../data/FinalData/Filtered_LSR_withNexgenCatchmentsgt1.csv")

plot_LSRincidents(ngen_cat, lsr_events)