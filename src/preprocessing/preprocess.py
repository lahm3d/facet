
# from whitebox_tools import WhiteboxTools
import whitebox
import geopandas as gpd
import pandas as pd
import fiona
from osgeo import gdal
from osgeo_utils.gdal_polygonize import gdal_polygonize
from pyproj.exceptions import ProjError
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
    preverving smoothing algorithm, and breaching any remaining depressions 
    (this step can be time consuming).

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
                filter = Config.preprocess['denoise']['filter_size'],
                norm_diff = Config.preprocess['denoise']['norm_diff'],
                num_iter = Config.preprocess['denoise']['num_iter'],
            )
            run_time = round( ( time.time() - start) / 60, 2 )
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

    # Breach depressions in the DEM using WhiteboxTools least cost depression breach algorithm
    if not Paths.breach.is_file():
        start = time.time()
        wbt.breach_depressions_least_cost(
            Paths.denoise,
            Paths.breach,
            dist=Config.preprocess['breach_depression_least_cost']['dist'],
            fill=Config.preprocess['breach_depression_least_cost']['fill'],
        )
        run_time = round( ( time.time() - start ) / 60, 2 )
        logger.info(f"Depressions breached. Run-time: {run_time} mins")
    else:
        logger.info("Depressions breached -- already exist!")


def create_weight_grid_from_streamlines(out_epsg, flowlines, watershed, dem,
                                        initiation_pixels, logger):
    """
    Creates a weight raster for TauDEM D8 flow accumulations using a pre-defined
    flowline file

    Parameters
    ----------
    out_epsg: Integer
        Integer specifying the European Petroleum Geospatial Group code defining
        the output horizontal coordinate reference system
    flowlines : WindowsPath object of pathlib module
        Path to the flowlines vector file (typically a *.shp).
    watershed : WindowsPath object of pathlib module
        Path to the watershed mask vector file (typically a *.shp).
    dem : WindowsPath object of pathlib module
        Path to the DEM raster file (typically a *.tif).
    initiation_pixels : WindowsPath object of pathlib module
        Path to the output stream initiation pixels raster file (typically a *.tif)..
    logger : Logger object of logging module
        Logger writes processing information to text file.

    Returns
    -------
    None.

    """
    if not initiation_pixels.is_file():
        # Read in the flowlines, clip them to the watershed mask
        flowlines = gpd.read_file(flowlines).to_crs( epsg = out_epsg)
        mask = utils.vector_to_geodataframe(watershed).to_crs( epsg = out_epsg )
        mask['geometry'] = mask.geometry.buffer( -1.0 ) # modify the watershed mask geometry to slightly smaller to prevent flowlines from having start nodes that are not inside the mask
        clip = gpd.clip(flowlines, mask)
        # multilinestrings get converted to linestrings
        clip = clip.explode( index_parts = True )

        # Generate end and start nodes for each linestring in the flowlines layer
        end_nodes = []
        start_nodes= []
        for line in clip['geometry']:
            end_nodes.append( tuple( np.array( line.coords )[-1][:2] ) )
            start_nodes.append( tuple( np.array( line.coords )[0][:2] ) )

        intersecting_nodes = set( start_nodes ).intersection( set( end_nodes ) ) # identify start nodes that are identical to end nodes
        
        # Create a GeoDataFrame from the start and end nodes lists
        init_nodes_geoms = []
        for node in start_nodes + end_nodes:
            if node not in intersecting_nodes:
                init_nodes_geoms.append( Point( node ) )
        init_nodes = gpd.GeoDataFrame( geometry = init_nodes_geoms, crs = flowlines.crs)

        # Open DEM to copy metadata and write a Weight Grid (WG) where 0 is not a channel initiation point and 1 is:
        with rasterio.open(dem) as src_dem:
            out_meta = src_dem.meta.copy()
            out_meta.update(compress = "lzw")
            out_meta.update(dtype = rasterio.int16)
            out_meta.update(nodata = -9999)

            with rasterio.open(initiation_pixels, "w+", **out_meta) as dst:
                array = dst.read(1)
                shapes = init_nodes['geometry'].values

                init_array = rasterio.features.rasterize(
                    shapes = shapes, default_value = 1, fill = 0, out = array,
                    transform = src_dem.transform,
                    # all_touched=True
                )

                init_array[init_array == -9999] = 0

                dst.write_band(1, init_array)
                logger.info("Channel initiation nodes generated")


def delineate_elevation_aligned_stream_network(Config, Paths, logger):
    """
    Run a series of TauDEM commands to generate a stream network
        1. d8flowdir - generates a D8 flow direction grid and a D8 slope grid 
        (drop/distance units) using the breach DEM as the input
        2. aread8 w/ weights - generates a D8 contributing area grid incorporating 
        weight grid representing channel initiation points from pre-defined flow
        network, checks for edge contamination
        3. aread8 w/o weights - generates a D8 contributing area grid without
        incorporating weight grid representing channel initiation points from 
        pre-defined flow network, checks for edge contamination
        4a. delineate watershed using weights - generates a vector stream network
        using the breach DEM, the D8 flow direction raster, weighted D8
        drainage area raster for drainage area, and the weighted D8 drainage area
        raster for the stream raster. Outputs a raster of Strahler order streams,
        a text file with a list links in the channel network tree, a textfile with a
        list of coordinates in the the channel network tree, an vector file of
        the channel network, and an output raster of watershed identifiers.
        4b. delineate watershed using a drainage area threshold - first, generate
        a stream raster by thresholding the unweighted D8 drainage area raster,
        then delineates the watershed similarly to 4A but using this new stream
        raster.
        5. dinf flow dir - generates d-infinity flow direction and slope rasters
        using the breach DEM.
        6a. dinf dist down using weights - generates a HAND grid using the non-hydro-enforced
        DEM, the d-infinity flow direction raster, and the weighted D8 drainage
        area raster
        6b. dinf dist down using threshold - generates a HAND grid using the non-hydro-enforced
        DEM, the d-infinity flow direction raster, and the stream raster generated
        using a drainage area threshold
        
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
    if not Paths.hand.is_file(): # Skip TauDEM processing if HAND outputs already exist in directory
        start = time.time()
        num_cores = Config.preprocess['taudem']['cores'] # set how many corse you want TauDEM to use
    
        # Create variables that will be fed to TauDEM commands using f-strings
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
        stream_raster = Paths.network_raster
        threshold = Config.preprocess['taudem']['drainage_threshold']
        
        taudem_workflow = {
            'd8 flow dir': f'mpiexec -n {num_cores} d8flowdir -fel "{breach}" -p "{p}" -sd8 "{sd8}"',
            'd8 flow acc w/ weights': f'mpiexec -n {num_cores} aread8 -p "{p}" -ad8 "{ad8_wg}" -wg "{init_pixels}" -nc',
            'd8 flow acc w/o weights': f'mpiexec -n {num_cores} aread8 -p "{p}" -ad8 "{ad8_no_wg}" -nc',
            'delineate watershed weights': (
                    f'''streamnet -fel "{breach}" -p "{p}" -ad8 "{ad8_no_wg}" -src "{ad8_wg}" -ord "{ord_g}" -tree "{tree}" -coord "{coord}" -net "{net}" -netlyr "{net.stem}" -w "{w}"'''
                ),
            'stream raster from d8': f'mpiexec -n {num_cores} threshold -ssa "{ad8_no_wg}" -src "{stream_raster}" -thresh "{threshold}" ',
            'delineate watershed threshold': (
            f'''streamnet -fel "{breach}" -p "{p}" -ad8 "{ad8_no_wg}" -src "{stream_raster}" -ord "{ord_g}" -tree "{tree}" -coord "{coord}" -net "{net}" -netlyr "{net.stem}" -w "{w}"'''
            ),
            'dinf flow dir': f'mpiexec -n {num_cores} dinfflowdir -fel "{breach}" -ang "{ang}" -slp "{slp}"',
            'hand grid weights': f'mpiexec -n {num_cores} dinfdistdown -fel "{dem}" -ang "{ang}" -src "{ad8_wg}" -dd "{dd}" -m ave v',
            'hand grid threshold': f'mpiexec -n {num_cores} dinfdistdown -fel "{dem}" -ang "{ang}" -src "{stream_raster}" -dd "{dd}" -m ave v'
            }
        
        # Assess which network delineation method is specified and modify TauDEM commands accordingly
        if Config.preprocess['taudem']['network_method'] == 'area_threshold':
            del taudem_workflow['delineate watershed weights']
            del taudem_workflow['hand grid weights']
        elif Config.preprocess['taudem']['network_method'] == 'flowline_weights':
            del taudem_workflow['stream raster from d8']
            del taudem_workflow['delineate watershed threshold']
            del taudem_workflow['hand grid threshold']
        else:
            logger.info("No stream network method specified!")
            raise ValueError("No stream network method specified!")
        # Run modified TauDEM commands        
        for cmd in taudem_workflow.values():
            run_command(cmd, logger)
        run_time = round( ( time.time() - start) / 60, 2 )
        logger.info( "TauDEM outputs created. Run-time: {run_time} mins" )            
    else:
        logger.info( "TauDEM outputs already exist. Skipping step." )

def polygonize_subwatershed_and_append_attributes(out_epsg, Config, Paths, logger):
    """
    Polygonizes the subwatersheds raster generated by TauDEM and appends
    attributes including physiography and information from the TauDEM-derived
    network.

    Parameters
    ----------
    out_epsg: Integer
        Integer specifying the European Petroleum Geospatial Group code defining
        the output horizontal coordinate reference system
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

    if not Paths.sub_watersheds_poly.is_file():
        # 1 - Polygonize the subwatersheds raster
        gdal_polygonize(
            str(Paths.sub_watersheds_rast),
            1,
            str(Paths.sub_watersheds_poly),
            connectedness8=True
        )

        # 2 - read in watershed polygons, use try/except block to deal with issues if local TauDEM install is missing some projection information (can happen on gov-controlled machines)
        try: 
            sub_watersheds = gpd.read_file(Paths.sub_watersheds_poly).to_crs( epsg = out_epsg )
        except ProjError:
            sub_watersheds = gpd.read_file(Paths.sub_watersheds_poly).set_crs( epsg = out_epsg, allow_override = True )
        sub_watersheds = sub_watersheds.rename( columns = {'DN': 'LINKNO'} )

        # 3 - polygons to centroid points
        watershed_points = sub_watersheds.copy()
        watershed_points.geometry = watershed_points.geometry.centroid

        # 4 - spatial join points to physiographic regions based on centroids
        physiography = utils.vector_to_geodataframe(
            Config.ancillary['physiography']
        ).to_crs( epsg = out_epsg )
        physio_in_watersheds = gpd.sjoin(watershed_points, physiography, how = 'left')  # spatial join
        physio_in_watersheds = physio_in_watersheds[ ['LINKNO', 'PROVINCE'] ]

        # 5 - merge physiography and network information to watershed polygons, use try/except block to deal with issues if local TauDEM install is missing some projection information (can happen on gov-controlled machines)
        try: 
            network = gpd.read_file(Paths.network_poly).to_crs( epsg = out_epsg )
        except (ProjError, ValueError):
            network = gpd.read_file(Paths.network_poly).set_crs( epsg = out_epsg, allow_override = True )
        network = network.drop( ['geometry'], axis = 1 )
        merge = sub_watersheds.merge( physio_in_watersheds, on = 'LINKNO' )  # merge 1 - physiography
        merge = merge.merge( network, on = 'LINKNO' )  # merge 2 - network information
        merge.to_file( Paths.sub_watersheds_poly )
        logger.info( "Sub-watersheds/catchments delineated" )
    else:
        logger.info( "Sub-watersheds/catchments already exist. Skipping step." )


def run_command(cmd, logger):
    """
    Execute commands as subprocesses

    Args:
        cmd: Command to run as a string
        logger: Logger instance

    Returns: None
    """
    try:
        p = subprocess.Popen( cmd, shell = True, stdout = subprocess.PIPE )
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
                a. merge and clip the crossings vetors (rails/roads) to the watershed polygon (mask)
                b. breach streams at crossings or burn cutlines
                c. denoise
                d. breach depressions
            3. create weight grid from streamlines
            4. TauDEM processing
                a. D8 flow direction
                b. D8 area accumulation 
                c. D8 area accumulation weighted by channel initiation points
                d. Channel network delineation (Channel initation weights or drainage area based)
                e. D-infinity flow direction
                h. HAND grid
            5. Vectorize subwatershed raster and attribute
    

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
    clip_flowlines( Config.ancillary['flowlines'], Paths.watershed, Paths.flowlines, logger )
    hydro_condition_dem( Config, Paths, logger )
    create_weight_grid_from_streamlines( int( Config.spatial_ref['epsg'] ), Paths.flowlines, Paths.watershed, Paths.dem, Paths.initiation_pixels, logger )
    delineate_elevation_aligned_stream_network( Config, Paths, logger )
    polygonize_subwatershed_and_append_attributes( int( Config.spatial_ref['epsg'] ), Config, Paths, logger )
