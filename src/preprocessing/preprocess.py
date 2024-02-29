
# from whitebox_tools import WhiteboxTools
import whitebox
import geopandas as gpd
import pandas as pd
import fiona
from osgeo import gdal
from osgeo_utils.gdal_polygonize import gdal_polygonize
import rasterio
from shapely.geometry import Point
import subprocess
import numpy as np
from utils import utils

def clip_flowlines(flowlines, mask, output, logger):

    aoi_flowlines = utils.vector_to_geodataframe(flowlines)
    mask = utils.vector_to_geodataframe(mask)
    flowlines = aoi_flowlines.clip(mask)

    try:
        flowlines.to_file(output)
        logger.info("NHD Flowlines clipped")
    except fiona.errors.DriverSupportError as e:
        logger.info(f"Error encountered while writing the flowline file: {e}")
        # datetime64 not supported by Esri Shapefile, so column will be dropped
        new_columns = [
            col for col in flowlines.columns if flowlines[col].dtype != 'datetime64[ms, UTC]'
            ]
        flowlines = flowlines[new_columns]
        flowlines.to_file(output)
        logger.info("NHD Flowlines clipped")


def merge_rails_and_roads(aoi_rails, aoi_roads, mask, output, logger):

    if not output.is_file():
        roads = utils.vector_to_geodataframe(aoi_roads)
        rails = utils.vector_to_geodataframe(aoi_rails)
        mask = utils.vector_to_geodataframe(mask)

        roads_mask = roads.clip(mask)
        rails_mask = rails.clip(mask)

        road_rail_crossings = gpd.GeoDataFrame(
            pd.concat([roads_mask, rails_mask], ignore_index=True), 
            crs=mask.crs
        )
        road_rail_crossings.to_file(output)
        logger.info("Roads and rails merged")

    logger.info("Roads and rails layer already exists. Skipping step")


def hydro_condition_dem(Config, Paths, logger):
    wbt = whitebox.WhiteboxTools()
    wbt._WhiteboxTools__compress_rasters = "True"
    wbt.set_verbose_mode(False)

    merge_rails_and_roads(
        Config.ancillary['census_rails'],
        Config.ancillary['census_roads'],
        Paths.watershed,
        Paths.road_rail_crossings,
        logger,
    )

    wbt.burn_streams_at_roads(
        Paths.dem, 
        Paths.flowlines, 
        Paths.road_rail_crossings, 
        Paths.burn_crossings, 
        width=Config.preprocess['burn_stream_at_roads']['width'], 
    )
    logger.info("Streams near roads burned")

    wbt.feature_preserving_smoothing(
        Paths.burn_crossings, 
        Paths.denoise, 
        filter=Config.preprocess['denoise']['filter_size'], 
        norm_diff=Config.preprocess['denoise']['norm_diff'], 
        num_iter=Config.preprocess['denoise']['num_iter'], 
    )
    logger.info("Feature preserving smoothing (denoising) performed")

    wbt.breach_depressions_least_cost(
        Paths.denoise,
        Paths.breach,
        dist=Config.preprocess['breach_depression_least_cost']['dist'],
        # max_cost=Config.preprocess['breach_depression_least_cost']['max_cost'],
        min_dist=None,
        fill=None,
    )


def create_weight_grid_from_streamlines(
    flowlines, watershed, dem, initiation_pixels
, logger):
    """
    Create weight file for TauDEM D8 FAC

    Args:
        flowlines:
        dem:
        initiation_pixels:

    Returns:
    """
    if not initiation_pixels.is_file():
        flowlines = gpd.read_file(flowlines)
        mask = utils.vector_to_geodataframe(watershed)

        mask['geometry'] = mask.geometry.buffer(-0.2)
        clip = gpd.clip(flowlines, mask)

        end_nodes = []
        start_nodes= []

        for line in clip['geometry']:
            end_nodes.append(tuple(np.array(line.coords)[-1][:2]))
            start_nodes.append(tuple(np.array(line.coords)[0][:2]))

        intersecting_nodes = set(start_nodes).intersection(set(end_nodes))

        init_nodes_geoms = []
        for node in start_nodes + end_nodes:
            if node not in intersecting_nodes:
                init_nodes_geoms.append(Point(node))
        init_nodes = gpd.GeoDataFrame(geometry=init_nodes_geoms, crs=flowlines.crs)

        # Open DEM to copy metadata and write a Weight Grid (WG):
        with rasterio.open(dem) as src_dem:
            out_meta = src_dem.meta.copy()
            out_meta.update(compress="lzw")
            out_meta.update(dtype=rasterio.int16)
            out_meta.update(nodata=-9999)

            with rasterio.open(initiation_pixels, "w+", **out_meta) as dst:
                array = dst.read(1)
                shapes = init_nodes['geometry'].values

                init_array = rasterio.features.rasterize(
                    shapes=shapes, default_value=1, fill=0, out=array, transform=src_dem.transform,
                    # all_touched=True
                )

                init_array[init_array == -9999] = 0

                dst.write_band(1, init_array)
                logger.info("Channel initiation nodes generated")



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

    # run_command(f'gdal_polygonize.py -8 "{Paths.sub_watersheds_rast}" "{Paths.sub_watersheds_poly}"')
    gdal_polygonize(
        str(Paths.sub_watersheds_rast),
        1,
        str(Paths.sub_watersheds_poly),
        connectedness8=True
    )

    # 1 - read in watershed polygons
    sub_watersheds = gpd.read_file(Paths.sub_watersheds_poly)
    sub_watersheds = sub_watersheds.rename(columns={'DN': 'LINKNO'})
    # 2 - polygons to centroid points
    watershed_points = sub_watersheds.copy()
    watershed_points.geometry = watershed_points.geometry.centroid
    # 3 - spatial join points to physiographic regions

    with fsspec.open(Config.ancillary['physiography']) as parquet:
        physiography = gpd.read_parquet(parquet)

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




def run_preprocessing_steps(Config, Paths, logger):
    clip_flowlines(Config.ancillary['flowlines'], Paths.watershed, Paths.flowlines, logger)
    hydro_condition_dem(Config, Paths, logger)
    create_weight_grid_from_streamlines(Paths.flowlines, Paths.watershed, Paths.dem, Paths.initiation_pixels, logger)
