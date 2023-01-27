



def setup_logging(hucID, huc_dir):
    # clean logging handlers
    clear_out_logger()

    # construct time stamps and log file name
    st_tstamp = strftime("%a, %d %b %Y %I:%M:%S %p") # e.g. Thu, 19 Sep 2019 12:20:41 PM
    log_name = hucID + strftime('_%y%m%d.log')
    log_file = huc_dir / log_name

    # delete old log file - it will not delete old logs from previous dates
    if log_file.is_file():
        os.remove(log_file)

    return log_file, st_tstamp


def initialize_logger(log_file):
    logger = logging.getLogger('logger_loader')
    logging.basicConfig(filename=log_file, filemode='a')
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        '%(asctime)s %(levelname)s [%(lineno)d] - %(message)s',
        '%m/%d/%Y %I:%M:%S %p'
        )
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


def clear_out_logger():
    # Remove all handlers associated with the root logger object.
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)


def open_memory_tif(arr, meta):
    from rasterio.io import MemoryFile
    #     with rasterio.Env(GDAL_CACHEMAX=256, GDAL_NUM_THREADS='ALL_CPUS'):
    with MemoryFile() as memfile:
        with memfile.open(**meta) as dataset:
            dataset.write(arr, indexes=1)
        return memfile.open()


def my_callback(value):
    if not "%" in value:
        print(value)


def get_cell_size(str_grid_path):
    with rasterio.open(str(str_grid_path)) as ds_grid:
        cs_x, cs_y = ds_grid.res
    return cs_x


def clean_tmp_files(huc_dir):
    """ cleans up intermediate files """

    # list of all files in HUC dir
    all_files = list(huc_dir.rglob('*'))

    keep_files = (
        list(huc_dir.rglob('*.log')) +
        list(huc_dir.rglob('*_dem.tif')) +
        list(huc_dir.rglob('*_mask.*')) +
        list(huc_dir.rglob('*_network.*')) +
        list(huc_dir.rglob('*_bankpixels.tif')) +
        list(huc_dir.rglob('*_floodplain.tif')) +
        list(huc_dir.rglob('*_hand.tif')) +
        list(huc_dir.rglob('*_channel_xns.*')) +
        list(huc_dir.rglob('post_processing/bankpoints_1D_metrics*.*')) +
        list(huc_dir.rglob('post_processing/floodplain_xns_1D_metrics*.*')) +
        list(huc_dir.rglob('post_processing/channel_floodplain_2D_metrics*.*'))
        )

    for file in all_files:
        if all([file.is_file(), file not in keep_files]):
            file.unlink()


def archive_huc(huc_dir):
    """ archives huc directory """

    # zip path
    huc_zip = Path(f'{huc_dir}.zip')

    if huc_zip.is_file():
        os.remove(huc_zip)
        shutil.make_archive(huc_dir, 'zip', huc_dir.parent, huc_dir.parts[-1])
    else:
        shutil.make_archive(huc_dir, 'zip', huc_dir.parent, huc_dir.parts[-1])


def reproject_ancillary_data(HUC04, PARAMS, logger):
    """
    Group re-projection of all ancillary datasets:
        NHD Streams
        NHD Waterbodies
        Physio file
        Census Roads
        Census Rails
    """

    # update PHYSIO param based on HUC04
    PHYSIO = PARAMS['physio drb'] if HUC04 == '0204' else PARAMS['physio cbw']

    reproject_list = [
        ('streams prj', PARAMS['ancillary dir'] / f'{HUC04}.shp'),
        ('waterbody prj', PARAMS['ancillary dir'] / f'{HUC04}_waterbody.shp'),
        ('physio prj', PHYSIO),
        ('census roads prj', PARAMS['census roads']),
        ('census rails prj', PARAMS['census rails']),
    ]

    for i in reproject_list:
        name, shp = i[0], Path(i[1])
        # check if native SHPs exist
        if shp.is_file():
            logger.info(f'\n{shp} file exists.')
            pass
        else:
            logger.error(f'\n{shp} file DOES NOT exist.')
            sys.exit(1)

        # reproject
        shp_proj = reproject_vector_layer(shp, PARAMS['crs'], logger)
        PARAMS[name] = Path(shp_proj)


# ==========================================================================
#   Reproject a grid layer using rasterio
# ==========================================================================
def define_grid_projection(str_source_grid, dst_crs, dst_file):
    print('Defining grid projection:')
    with rasterio.open(str_source_grid, 'r') as src:
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
        })
        arr_src = src.read(1)

    with rasterio.open(dst_file, 'w', **kwargs) as dst:
        dst.write(arr_src, indexes=1)


# ==========================================================================
#   Reproject a grid layer using rasterio
# ==========================================================================
def reproject_grid_layer(str_source_grid, dst_crs, dst_file, resolution, logger):
    # reproject raster plus resample if needed
    # Resolution is a pixel value as a tuple
    try:
        st = timer()
        with rasterio.open(str_source_grid) as src:
            transform, width, height = calculate_default_transform(
                src.crs, dst_crs, src.width, src.height, *src.bounds, resolution=resolution)
            kwargs = src.meta.copy()
            kwargs.update({
                'crs': dst_crs,
                'transform': transform,
                'width': width,
                'height': height
            })

            with rasterio.open(dst_file, 'w', **kwargs) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, i),
                        destination=rasterio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=dst_crs,
                        resampling=Resampling.bilinear)
        end = round((timer() - st)/60.0, 2)
        logger.info(f'Reprojected DEM. Time elapsed: {end} mins')
        return dst_file
    except:
        logger.critical(f'{str_source_grid}: failed to reproject.')
        sys.exit(1)


# ==========================================================================
#   Reproject a vector layer using geopandas
# ==========================================================================
def reproject_vector_layer(in_shp, str_target_proj4, logger):
    print(f'Reprojecting vector layer: {in_shp}')

    proj_shp = in_shp.parent / f'{in_shp.stem}_proj.shp'

    if proj_shp.is_file():
        logger.info(f'{proj_shp} reprojected file already exists\n')
        return str(proj_shp)
    else:
        gdf = gpd.read_file(str(in_shp))
        # fix float64 to int64
        float64_2_int64 = ['NHDPlusID', 'Shape_Area', 'DSContArea', 'USContArea']
        for col in float64_2_int64:
            try:
                gdf[col] = gdf[col].astype(np.int64)
            except KeyError:
                pass

        gdf_proj = gdf.to_crs(str_target_proj4)
        gdf_proj.to_file(str(proj_shp))

        logger.info(f'{proj_shp} successfully reprojected\n')
        return str(proj_shp)


# ==========================================================================
#   For clipping features
# ==========================================================================
def clip_features_using_grid(
        str_lines_path, output_filename, str_dem_path, in_crs, logger, mask_shp):
    # clip features using HUC mask, if the mask doesn't exist polygonize DEM
    mask_shp = Path(mask_shp)
    if mask_shp.is_file():
        st = timer()
        # whitebox clip
        WBT.clip(str_lines_path, mask_shp, output_filename)
        end = round((timer() - st)/60.0, 2)
        logger.info(f'Streams clipped by {mask_shp}. Time Elapsed: {end} mins')
    else:
        st = timer()
        logger.warning(f'''
        {mask_shp} does not file exists. Creating new mask from DEM.
        This step can be error prone please review the output.
        ''')
        # Polygonize the raster DEM with rasterio:
        with rasterio.open(str(str_dem_path)) as ds_dem:
            arr_dem = ds_dem.read(1)

        arr_dem[arr_dem > 0] = 100
        mask = arr_dem == 100

        results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) in enumerate(shapes(arr_dem, mask=mask, transform=ds_dem.transform))
            )

        poly = list(results)
        poly_df = gpd.GeoDataFrame.from_features(poly)
        poly_df.crs = in_crs
        # poly_df = poly_df[poly_df.raster_val == 100.0]
        # tmp_shp = os.path.dirname(str_dem_path) + "/mask.shp"  # tmp huc mask
        poly_df.to_file(str(mask_shp))

        # whitebox clip
        WBT.clip(str_lines_path, str(mask_shp), output_filename)
        end = round((timer() - st)/60.0, 2)
        logger.info(f'Streams clipped by {mask_shp}. Time Elapsed: {end} mins')


# ==========================================================================
#   Polygonize watersheds
# ==========================================================================
def watershed_polygonize(in_tif, out_shp, dst_crs, logger):
    logger.info("Polygonizing reach catchments")
    tmp = os.path.dirname(in_tif) + "\\breach_w_tmp.shp"  # tmp DEM mask

    # convert the raster to polygon
    mask = None
    with rasterio.open(in_tif) as src:
        image = src.read(1)
        results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) in enumerate(shapes(image, mask=mask, transform=src.transform))
            )

    # write raster vals to polygon
    driver = 'Shapefile'
    crs = dst_crs
    schema = {'properties': [('raster_val', 'int')], 'geometry': 'Polygon'}
    with fiona.open(str(tmp), 'w', driver=driver, crs=crs, schema=schema) as dst:
        dst.writerecords(results)

    # original sauce: https://gis.stackexchange.com/questions/149959/dissolving-polygons-based-on-attributes-with-python-shapely-fiona
    # dissolve and remove nodata vals
    with fiona.open(tmp) as input:
        meta = input.meta  # copy tmp files metadata
        with fiona.open(out_shp, 'w', **meta) as output:
            # sort and then perform groupby on raster_vals
            by_vals = sorted(input, key=lambda k: k['properties']['raster_val'])
            # group by 'raster_val'
            for key, group in itertools.groupby(
                    by_vals, key=lambda x: x['properties']['raster_val']):
                properties, geom = zip(*[(feature['properties'], shape(feature['geometry'])) for feature in group])
                # perform check and exclude nodata value
                if properties[0]['raster_val'] >= 0:
                    # capture & exclude geometries that throw ValueError performing geometry union
                    try:
                        geom = [g if g.is_valid else g.buffer(0.0) for g in geom]
                        test_unary = unary_union(geom)
                        # write the feature, computing the unary_union of the elements...
                        # ...in the group with the properties of the first element in the group
                        output.write(
                            {'geometry': mapping(unary_union(geom)), 'properties': properties[0]}
                            )
                    except ValueError:
                        logger.warning(f'''
                        catchment value ({properties[0]['raster_val']})
                        encountered a topology exception error.
                        Attempting to fix geometry...
                        ''')
                        """ if union fails then loop through geom,
                        simplify geometry, append to a new geom and then perform
                        union -- should resolve most self-intersecting points"""
                        new_geom = []
                        for g in geom:
                            g_2_list = mapping(g)['coordinates'][0]
                            simplified_geom_lst = list(Polygon(g_2_list).simplify(0).exterior.coords)
                            new_geom.append(Polygon(simplified_geom_lst))
                        # write out updated geometry
                        output.write(
                            {'geometry': mapping(unary_union(new_geom)), 'properties': properties[0]}
                            )
                        # log warning about fixed geometry
                        logger.warning(f'''
                        catchment value ({properties[0]['raster_val']}) topology error fixed.
                        Please review catchment file.
                        ''')


# ==========================================================================
#   join watersheds attributes
# ==========================================================================
def join_watershed_attrs(w, physio, net, output, logger):
    # 1 - read in watershed polygons
    wsheds = gpd.read_file(str(w))  # watershed grid shp
    points = wsheds.copy()  #
    # 2 - polygons to centroid points
    points.geometry = points['geometry'].centroid  # get centroids
    points.crs = wsheds.crs  # copy poly crs
    # 3 -  spatial join points to physiographic regions
    physio = gpd.read_file(str(physio))
    pointsInPhysio = sjoin(points, physio, how='left')  # spatial join

    # 4 - merge attrs to watershed polygons
    net = gpd.read_file(str(net))

    # 4.1 - rename columns
    wsheds = wsheds.rename(columns={'raster_val': 'LINKNO'})
    pointsInPhysio = pointsInPhysio.rename(columns={'raster_val': 'LINKNO'})

    # 4.2 - drop geometry from pointsInPhysio & net files
    pointsInPhysio = pointsInPhysio.drop(['geometry'], axis=1)
    net = net.drop(['geometry'], axis=1)

    # 4.3 - merge pointsInPhysio & net files to wsheds
    wshedsMerge = wsheds.merge(pointsInPhysio, on='LINKNO')  # merge 1
    wshedsMerge = wshedsMerge.merge(net, on='LINKNO')  # merge 2

    wshedsMerge.to_file(str(output))
    logger.info('Stream and Physio attributes successfully joined to Catchments')


def rasterize_gdf(str_net_path, str_ingrid_path, str_tempgrid):
    gdf = gpd.read_file(str(str_net_path))
    '''
    Thanks to:
    https://gis.stackexchange.com/questions/151339/rasterize-a-shapefile-with-geopandas-or-fiona-python
    '''
    with rasterio.open(str_ingrid_path) as rst:
        meta = rst.meta.copy()
    meta.update(compress='lzw')
    meta.update(dtype=rasterio.int32)
    meta.update(nodata=0)

    with rasterio.open(str_tempgrid, 'w+', **meta) as out:
        out_arr = out.read(1)
        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom, value) for geom, value in zip(gdf.geometry, gdf['LINKNO']))
        arr_burned = rasterio.features.rasterize(
            shapes=shapes, fill=0, out=out_arr, transform=out.transform
            )
        out.write_band(1, arr_burned)

    return


# Count the number of features in a vector file:
def get_feature_count(str_shp_path):
    with fiona.open(str(str_shp_path), 'r') as features:
        i_count = len(features)
    return i_count
