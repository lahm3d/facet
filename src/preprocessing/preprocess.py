
# from whitebox_tools import WhiteboxTools
import whitebox
import geopandas as gpd
import pandas as pd
import fiona
from osgeo import gdal
from osgeo_utils.gdal_polygonize import gdal_polygonize
import rasterio
import rasterstats
from shapely.geometry import Point
import subprocess
import numpy as np
from src.utils import utils
import time

def clip_flowlines(flowlines, mask, output, logger):
    """
    Clips the provided flowlines to the watershed polygon boundary in case they
    extend beyond it which would create problems in subsequent steps

    Parameters
    ----------
    flowlines : str
        String defining the filepath for the input flowlines *.shp file (typically NHD Plus HR flowlines).
    mask : WindowsPath object of pathlib module
        Path to watershed polygon (typically a *.shp file).
    output : WindowsPath object of pathlib module
        Path where clipped flowlines are written (typically a *.shp file).
    logger : Logger object of logging module
        Logger writes processing information to text file.

    Returns
    -------
    None.

    """

    if not output.is_file(): # only enter this step if the file does not already exist
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


def merge_rails_and_roads(out_epsg, aoi_rails, aoi_roads, mask, output, logger):
    """
    Clips the roads and rails line vector files to the watershed mask, and then
    merges the clipped line vectors to a single combined file which is written
    to the data/version directory.

    Parameters
    ----------
    out_epsg : Integer
        Integer specifying the European Petroleum Geospatial Group code defining
        the output horizontal coordinate reference system
    aoi_rails : String
        String defining the filepath for the railroad lines vector file (typically a *.shp file).
    aoi_roads : String
        String defining the filepath for the road lines vector file (typically a *.shp file).
    mask : WindowsPath object of pathlib module
        Path to watershed polygon (typically a *.shp file).
    output : WindowsPath object of pathlib module
        Path to watershed polygon (typically a *.shp file).
    logger : Logger object of logging module
        Logger writes processing information to text file.

    Returns
    -------
    None.

    """

    if not output.is_file(): # if the merged roads/rails file already exists, do not proceed with this step
        roads = utils.vector_to_geodataframe(aoi_roads)
        rails = utils.vector_to_geodataframe(aoi_rails)
        mask = utils.vector_to_geodataframe(mask).to_crs( epsg = out_epsg )

        roads_mask = roads.clip( mask.to_crs( roads.crs ) ).to_crs( epsg = out_epsg )
        rails_mask = rails.clip( mask.to_crs( rails.crs ) ).to_crs( epsg = out_epsg )

        road_rail_crossings = gpd.GeoDataFrame(
            pd.concat( [ roads_mask, rails_mask ], ignore_index = True), 
            crs = mask.crs
        )
        road_rail_crossings.to_file(output)
        logger.info("Roads and rails merged")
    else:
        logger.info("Roads and rails layer already exists. Skipping step")


def burn_cutlines(dem, output, cutlines, out_epsg, logger):
    """
    Burn the cutlines into the DEM

    Parameters
    ----------
    dem : WindowsPath object of pathlib module
        Path to the input DEM which will be cut.
    output : WindowsPath object of pathlib module
        Path to which the cut DEM will be written.
    cutlines : String
        Path to the cutlines vector file (typically a *.
    out_epsg : Integer
        DESCRIPTION.
    logger : Logger object of logging module
        Logger writes processing information to text file.

    Returns
    -------
    None.

    """
    
    # Read in the cutlines and attribute them with the minimum DEM elevation
    cutlines_gdf = utils.vector_to_geodataframe(cutlines).to_crs( epsg = out_epsg )
    cut_zmin = rasterstats.zonal_stats( cutlines, dem, stats = "min" )
    cutlines_gdf['elevation'] = [i['min'] for i in cut_zmin]
    cutline_zip = zip(cutlines_gdf.geometry.values, cutlines_gdf.elevation.values)
    
    # Burn the cutlines into the DEM and write out to a new file
    with rasterio.open(dem, "r") as src:
        band_num = 1
        src_image = src.read(band_num, out_dtype = src.meta['dtype'])
        dem_cut = rasterio.features.rasterize( cutline_zip, out = src_image,
                                              transform = src.transform,
                                              all_touched = True )
        
        # save tif
        profile = src.profile
        profile.update( dtype= src.meta['dtype'], count = 1, compress = "lzw" )

        with rasterio.open( output, "w", **profile ) as dst:
            dst.write( dem_cut, 1 )


def hydro_condition_dem(Config, Paths, logger):
    """
    Hydro-condition the provided DEM by burning in cutlines, burning in road/rail
    crossings near pre-defined streams, denoising the DEM by applying a feature-
    preverving smoothing algorithm, and breaching any remaining depressions.

    Parameters
    ----------
    Config : CreatConfig object of src.utils.parse_tomle module
        Object containing dictionaries specifying processing parameters.
    Paths : CreateFilepaths object of src.utils.parse_tomle module
        Object containing filepaths for the outputs of the FACET workflow.
    logger : Logger object of logging module
        Logger writes processing information to text file.
        
    Returns
    -------
    None.

    """
    
    # Setup whiteboxtool options
    wbt = whitebox.WhiteboxTools()
    wbt._WhiteboxTools__compress_rasters = "True"
    wbt.set_verbose_mode(False)

    # Merge rail/road features into a single file for "burn_stream_at_roads" function
    merge_rails_and_roads(
        int( Config.spatial_ref['epsg'] ),
        Config.ancillary['census_rails'],
        Config.ancillary['census_roads'],
        Paths.watershed,
        Paths.road_rail_crossings,
        logger,
    )
    
    # Burn cutlines into DEM
    if ( not Paths.burn_cutlines.is_file() ) & ( Config.preprocess['burn_cutlines']['burn_cutlines_flag'] == True ):
        start = time.time()
        burn_cutlines( Paths.dem, Paths.burn_cutlines, Config.ancillary['cutlines'], int( Config.spatial_ref['epsg'] ), logger )
        run_time = round((time.time() - start) / 60, 2)
        logger.info(f"Cutlines burned. Run-time: {run_time} mins")
    elif Config.preprocess['burn_cutlines']['burn_cutlines_flag'] == False:
        logger.info("Cutlines not burned -- burn_cutlines_flag == False")
    else:
        logger.info("Cutlines burned -- already exist!")


    # Burn rail/road crossings into DEM where they intersect the flowlines vector file
    if ( not Paths.burn_crossings.is_file() ) & ( Config.preprocess['burn_stream_at_roads']['burn_stream_at_roads_flag'] == True ):
        start = time.time()
        wbt.burn_streams_at_roads(
            Paths.dem, 
            Paths.flowlines, 
            Paths.road_rail_crossings, 
            Paths.burn_crossings, 
            width = Config.preprocess['burn_stream_at_roads']['width'], 
        )
        run_time = round((time.time() - start) / 60, 2)
        logger.info(f"Streams near roads burned. Run-time: {run_time} mins")
    elif  Config.preprocess['burn_stream_at_roads']['burn_stream_at_roads_flag'] == False:
        logger.info("Streams near roads not burned -- burn_stream_at_roads_flag == False!")
    else:
        logger.info("Streams near roads burned -- already exist!")

    # Denoise the DEM, or change the Paths.denoise path if the denoise_flag == False
    if Config.preprocess['denoise']['denoise_flag'] == True: # Only enter the denoise process if the denoise_flag is True
        if not Paths.denoise.is_file():
            start = time.time()
            wbt.feature_preserving_smoothing(
                Paths.burn_crossings, 
                Paths.denoise, 
                filter=Config.preprocess['denoise']['filter_size'],
                norm_diff=Config.preprocess['denoise']['norm_diff'],
                num_iter=Config.preprocess['denoise']['num_iter'],
            )
            run_time = round((time.time() - start) / 60, 2)
            logger.info(f"Feature preserving smoothing (denoising) performed. Run-time: {run_time} mins")
        else:
            logger.info("Feature preserving smoothing (denoising) performed -- already exist!")
    else:
        if Paths.burn_cutlines.is_file():
            Paths.denoise = Paths.burn_cutlines
            logger.info("Feature preserving smoothing (denoising) not performed -- flag = False, using the unsmoothed, cutline burned DEM for future steps!")
        if Paths.burn_crossings.is_file():
            Paths.denoise = Paths.burn_crossings
            logger.info("Feature preserving smoothing (denoising) not performed -- flag = False, using the unsmoothed, rail/road crossings burned DEM for future steps!")
     
    if not Paths.breach.is_file():
        start = time.time()
        wbt.breach_depressions_least_cost(
            Paths.denoise,
            Paths.breach,
            dist=Config.preprocess['breach_depression_least_cost']['dist'],
            fill=Config.preprocess['breach_depression_least_cost']['fill'],
        )
        run_time = round((time.time() - start) / 60, 2)
        logger.info(f"Depressions breached. Run-time: {run_time} mins")
    else:
        logger.info("Depressions breached -- already exist!")


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

        mask['geometry'] = mask.geometry.buffer(-1.0)
        clip = gpd.clip(flowlines, mask)
        # multilinestrings get converted to linestrings
        clip = clip.explode(index_parts=True)

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


def delineate_elevation_aligned_stream_network(Config, Paths, logger):

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

    for cmd in taudem_workflow.values():
        run_command(cmd, logger)
    
    logger.info("Preprocessing steps complete")


def polygonize_subwatershed_and_append_attributes(Config, Paths, logger):

    if not Paths.sub_watersheds_poly.is_file():
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
        physiography = utils.vector_to_geodataframe(
            Config.ancillary['physiography']
        )

        physio_in_watersheds = gpd.sjoin(watershed_points, physiography, how='left')  # spatial join
        physio_in_watersheds = physio_in_watersheds[['LINKNO', 'PROVINCE']]

        # 4 - merge attrs to watershed polygons
        network = gpd.read_file(Paths.network_poly)
        network = network.drop(['geometry'], axis=1)

        merge = sub_watersheds.merge(physio_in_watersheds, on='LINKNO')  # merge 1
        merge = merge.merge(network, on='LINKNO')  # merge 2
        merge.to_file(Paths.sub_watersheds_poly)
        logger.info("Sub-watersheds/catchments delineated")
    else:
        logger.info("Sub-watersheds/catchments already exist. Skipping step.")


def run_command(cmd, logger):
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
            logger.error(f'\n {text} \n')
        else:
            logger.info(f'Run successful. \n {text}')

    except subprocess.CalledProcessError as e:
        logger.error(f"failed to return code: {e}")
    except OSError as e:
        logger.error(f"failed to execute shell: {e}")
    except IOError as e:
        logger.error(f"failed to read file(s): {e}")


def run_preprocessing_steps(Config, Paths, logger):
    """ 
    Run the pre-processing steps:
            1. clip flowlines to watershed polygon
            2. hydro-condition the DEM
                a.
                b.
                c.
            3. create weight grid from streamlines
            4.
    

    Parameters
    ----------
    Config : CreatConfig object of src.utils.parse_tomle module
        Object containing dictionaries specifying processing parameters.
    Paths : CreateFilepaths object of src.utils.parse_tomle module
        Object containing filepaths for the outputs of the FACET workflow.
    logger : Logger object of logging module
        Logger writes processing information to text file.

    Returns
    -------
    None.

    """
    clip_flowlines(Config.ancillary['flowlines'], Paths.watershed, Paths.flowlines, logger)
    hydro_condition_dem(Config, Paths, logger)
    create_weight_grid_from_streamlines(Paths.flowlines, Paths.watershed, Paths.dem, Paths.initiation_pixels, logger)
    delineate_elevation_aligned_stream_network(Config, Paths, logger)
    polygonize_subwatershed_and_append_attributes(Config, Paths, logger)