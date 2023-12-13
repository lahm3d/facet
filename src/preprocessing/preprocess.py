
# from whitebox_tools import WhiteboxTools
import whitebox
import geopandas as gpd
import pandas as pd
import fiona
import rasterio
from rasterio.warp import transform
from rasterio.features import shapes
from shapely.geometry import shape
import subprocess
import numpy as np

def clip_flowlines(flowlines, mask, output):

    aoi_flowlines = gpd.read_parquet(flowlines)
    mask = gpd.read_file(mask)
    flowlines = aoi_flowlines.clip(mask)

    # to address datetime64 not supported by Esri Shapefile
    datetime_columns = [
        col for col in flowlines.columns if flowlines[col].dtype == 'datetime64[ms, UTC]'
        ]
    for col in datetime_columns:
        flowlines[col] = flowlines[col].dt.strftime('%Y-%m-%d %H:%M:%S')

    flowlines.to_file(output)


def merge_rails_and_roads(aoi_rails, aoi_roads, mask, output):

    roads = gpd.read_parquet(aoi_roads)
    rails = gpd.read_parquet(aoi_rails)

    mask = gpd.read_file(mask)

    roads_mask = roads.clip(mask)
    rails_mask = rails.clip(mask)

    road_rail_crossings = gpd.GeoDataFrame(
        pd.concat([roads_mask, rails_mask], ignore_index=True), 
        crs=mask.crs
    )
    road_rail_crossings.to_file(output)


def hydro_condition_dem(Config, Paths):
    wbt = whitebox.WhiteboxTools()
    wbt._WhiteboxTools__compress_rasters = "True"
    wbt.set_verbose_mode(False)

    merge_rails_and_roads(
        Config.ancillary['census_rails'], 
        Config.ancillary['census_roads'], 
        Paths.watershed, 
        Paths.road_rail_crossings
    )

    wbt.burn_streams_at_roads(
        Paths.dem, 
        Paths.flowlines, 
        Paths.road_rail_crossings, 
        Paths.burn_crossings, 
        width=Config.preprocess['burn_stream_at_roads']['width'], 
    )

    wbt.feature_preserving_smoothing(
        Paths.burn_crossings, 
        Paths.denoise, 
        filter=Config.preprocess['denoise']['filter_size'], 
        norm_diff=Config.preprocess['denoise']['norm_diff'], 
        num_iter=Config.preprocess['denoise']['num_iter'], 
    )

    wbt.breach_depressions_least_cost(
        Paths.denoise, 
        Paths.breach, 
        dist=Config.preprocess['breach_depression_least_cost']['dist'], 
        fill=Config.preprocess['breach_depression_least_cost']['fill'], 
    )


def create_weight_grid_from_streamlines(
    flowlines, dem, initiation_pixels
):
    """
    Create weight file for TauDEM D8 FAC

    Args:
        flowlines:
        dem:
        initiation_pixels:

    Returns:
    """

    lst_coords = []
    lst_pts = []
    lst_x = []
    lst_y = []

    with fiona.open(flowlines) as lines:

        streamlines_crs = lines.crs  # to use in the output grid

        # Get separate lists of start and end points:
        for line in lines:
            if line["geometry"]["type"] == "LineString":  # Make sure it's a LineString
                # Add endpts:
                lst_coords.append(line["geometry"]["coordinates"][-1])
                # Add startpts:
                lst_pts.append(line["geometry"]["coordinates"][0])

        # If a start point is not also in the endpt list, it's first order:
        for pt in lst_pts:
            if pt not in lst_coords:
                lst_x.append(pt[0])
                lst_y.append(pt[1])

    # Open DEM to copy metadata and write a Weight Grid (WG):
    with rasterio.open(dem) as ds_dem:
        out_meta = ds_dem.meta.copy()
        out_meta.update(compress="lzw")
        out_meta.update(dtype=rasterio.int16)
        out_meta.update(nodata=-9999)
        # out_meta.update(crs=lines.crs)  # shouldn't be necessary

        # Construct the output array:
        arr_danglepts = np.zeros(
            [out_meta["height"], out_meta["width"]], dtype=out_meta["dtype"]
        )

        tpl_pts = transform(streamlines_crs, out_meta["crs"], lst_x, lst_y)
        lst_dangles = zip(tpl_pts[0], tpl_pts[1])

        for coords in lst_dangles:
            # BUT you have to convert coordinates from hires to dem
            col, row = ~ds_dem.transform * (coords[0], coords[1])
            try:
                arr_danglepts[int(row), int(col)] = 1
            except:
                continue

    # Now write the new grid using this metadata:
    with rasterio.open(initiation_pixels, "w", **out_meta) as dest:
        dest.write(arr_danglepts, indexes=1)


def delineate_elevation_aligned_stream_network(Config, Paths):

    num_cores = Config.preprocess['taudem']['cores']

    breach = Paths.breach
    p = Paths.d8_fdir_point
    sd8 = Paths.slope_grid_sd8
    ad8_wg = Paths.area_grid_ad8_ip
    init_pixels = Paths.initiation_pixels
    ad8_no_wg = Paths.area_grid_ad8
    ord_g = Paths.network_order
    tree = Paths.network_tree
    coord = Paths.network_coords
    net = Paths.network_poly
    w = Paths.sub_watersheds_rast
    ang = Paths.flow_dir_dinf
    slp = Paths.slope_grid_dinf
    dem = Paths.dem
    dd = Paths.hand
    
    taudem_workflow = {
        'd8 flow dir': f'mpiexec -n {num_cores} d8flowdir -fel "{breach}" -p "{p}" -sd8 "{sd8}"',
        'd8 flow acc w/ grid': f'mpiexec -n {num_cores} aread8 -p "{p}" -ad8 "{ad8_wg}" -wg "{init_pixels}" -nc',
        'd8 flow acc w/o grid': f'mpiexec -n {num_cores} aread8 -p "{p}" -ad8 "{ad8_no_wg}" -nc',
        'delineate watershed': (
                f'''streamnet -fel "{breach}" -p "{p}" -ad8 "{ad8_no_wg}" -src "{ad8_wg}" -ord "{ord_g}" -tree "{tree}" -coord "{coord}" -net "{net}" -netlyr "{net.stem}" -w "{w}"'''
            ),
        'dinf flow dir': f'mpiexec -n {num_cores} dinfflowdir -fel "{breach}" -ang "{ang}" -slp "{slp}"',
        'hand grid': f'mpiexec -n {num_cores} dinfdistdown -fel "{dem}" -ang "{ang}" -src "{ad8_wg}" -dd "{dd}" -m ave v'
        }

    for key, cmd in taudem_workflow.items():
        run_command(cmd)

def polygonize_subwatershed_and_append_attributes(Config, Paths):

    run_command(f'gdal_polygonize.py -8 "{Paths.sub_watersheds_rast}" "{Paths.sub_watersheds_poly}"')

    # 1 - read in watershed polygons
    sub_watersheds = gpd.read_file(Paths.sub_watersheds_poly)
    sub_watersheds = sub_watersheds.rename(columns={'DN': 'LINKNO'})
    # 2 - polygons to centroid points
    watershed_points = sub_watersheds.copy()
    watershed_points.geometry = watershed_points.geometry.centroid
    # 3 - spatial join points to physiographic regions
    physiography = gpd.read_parquet(Config.ancillary['physiography'])
    physio_in_watersheds = gpd.sjoin(watershed_points, physiography, how='left')  # spatial join
    physio_in_watersheds = physio_in_watersheds[['LINKNO', 'PROVINCE']]

    # 4 - merge attrs to watershed polygons
    network = gpd.read_file(Paths.network_poly)
    network = network.drop(['geometry'], axis=1)

    merge = sub_watersheds.merge(physio_in_watersheds, on='LINKNO')  # merge 1
    merge = merge.merge(network, on='LINKNO')  # merge 2
    merge.to_file(Paths.sub_watersheds_poly)


def run_command(cmd: str) -> None:
    """
    Execute commands as subprocesses

    Args:
        cmd: Command to run as a string
        logger: Logger instance

    Returns: None
    """
    try:
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        output, err = p.communicate()

        # Get some feedback from the process to print out:
        if err is None:
            text = output.decode()
            print("\n", text, "\n")
        else:
            print(err)

    except subprocess.CalledProcessError as e:
        print(f"failed to return code: {e}")
    except OSError as e:
        print(f"failed to execute shell: {e}")
    except IOError as e:
        print(f"failed to read file(s): {e}")


def run_preprocessing_steps(Config, Paths):
    clip_flowlines(Config.ancillary['flowlines'], Paths.watershed, Paths.flowlines)
    hydro_condition_dem(Config, Paths)
    create_weight_grid_from_streamlines(Paths.flowlines, Paths.dem, Paths.initiation_pixels)
    delineate_elevation_aligned_stream_network(Config, Paths)
    polygonize_subwatershed_and_append_attributes(Config, Paths)

# if __name__ == '__main__':




