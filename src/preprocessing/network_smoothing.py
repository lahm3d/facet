import numpy as np
import geopandas as gpd
from shapely.geometry import LineString

def chaikins_corner_cutting(geom, refinements=5):
    coords = np.array(geom.coords)

    for _ in range(refinements):
        L = coords.repeat(2, axis=0)
        R = np.empty_like(L)
        R[0] = L[0]
        R[2::2] = L[1:-1:2]
        R[1:-1:2] = L[2::2]
        R[-1] = L[-1]
        coords = L * 0.75 + R * 0.25
    
    return LineString(coords)

def calculate_straight_line_length(geom):
    start_node = geom.interpolate(0)
    end_node = geom.interpolate(1, normalized=True)
    return start_node.distance(end_node)

def apply_chaikins_corner_cutting(input_file, output_file, refinements=5):
    gdf = gpd.read_file(input_file)
    gdf['geometry'] = gdf['geometry'].apply(lambda geom: chaikins_corner_cutting(geom, refinements))
    gdf['sm_length'] = gdf['geometry'].length
    gdf['sm_SL'] = gdf['geometry'].apply(lambda geom: calculate_straight_line_length(geom))
    gdf.to_file(output_file)