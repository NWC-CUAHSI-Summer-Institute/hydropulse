import geopandas as gpd
from pathlib import Path

def calculate_area(catchment_id, catchments_dir):
    """Calculate area from a single GeoPackage file for the given catchment."""
    gpkg_file = Path(catchments_dir) / "nex_catchments_subset.gpkg"
    if not gpkg_file.exists():
        raise FileNotFoundError(f"GeoPackage file not found at {gpkg_file}")
    
    # Load the 'divides' layer from the GeoPackage
    gdf = gpd.read_file(gpkg_file, layer="divides")
    
    # Filter for the specific catchment_id
    catchment_gdf = gdf[gdf["divide_id"] == catchment_id]
    if catchment_gdf.empty:
        raise ValueError(f"No geometry found for catchment_id {catchment_id}")
    
    # Reproject to a projected CRS (in meters) if not already projected
    if not catchment_gdf.crs.is_projected:
        catchment_gdf = catchment_gdf.to_crs(epsg=3857)  # Web Mercator
    
    # Calculate area in square kilometers
    area_sqkm = catchment_gdf.geometry.area.iloc[0] / 1e6
    return area_sqkm