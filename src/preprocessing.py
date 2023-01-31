def breach_dem(str_dem_path, breach_output, logger):
    # breach dem
    st = timer()
    WBT.breach_depressions(str_dem_path, breach_output, fill_pits=False)
    end = round((timer() - st) / 60.0, 2)
    logger.info(
        f"{breach_output}: breach output successfully created. Time elapsed: {end} mins"
    )


def run_tauDEM(cmd, logger):
    """execute tauDEM commands"""

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
        logger.critical(f"failed to return code: {e}")
    except OSError as e:
        logger.critical(f"failed to execute shell: {e}")
    except IOError as e:
        logger.critical(f"failed to read file(s): {e}")


# ===============================================================================
#  Hydrologically condition DEM to allow breaching road & stream x-cross sections
# ===============================================================================
def pre_breach_DEM_conditioning(
    huc_dir, hucID, str_dem_path, str_nhd_path, census_roads, census_rails, mask, logger
):
    # pre-breach DEM conditioning

    st = timer()
    # construct paths for temp & final files
    tmp_roads = huc_dir / "tmp_roads.shp"
    tmp_rails = huc_dir / "tmp_rails.shp"
    x_sect_pts = str(huc_dir / "x_section_pts.shp")
    x_sect_polys = str(huc_dir / "x_section_polys.shp")
    ds_min_filter = str(huc_dir / "ds_min_filter.tif")
    ds_min_clip = str(huc_dir / "ds_min_clip.tif")
    dem_merge = str(huc_dir / "dem_road_stream.tif")

    # clip census roads and rails by HUC mask
    shps = [(census_roads, str(tmp_roads)), (census_rails, str(tmp_rails))]

    for shp in shps:
        # whitebox clip
        inSHP, outSHP = shp[0], shp[1]
        WBT.clip(inSHP, mask, outSHP, callback=my_callback)

    # Read line layers:
    gdf_nhd = gpd.read_file(str(str_nhd_path))
    gdf_roads = gpd.read_file(str(tmp_roads))
    if tmp_rails.is_file():
        gdf_rails = gpd.read_file(str(tmp_rails))

        ###---If rail roads exist within the HUC---###

        # find unary_unions as points for streams x roads and rails
        pts_0 = gdf_roads.unary_union.intersection(gdf_nhd.unary_union)
        pts_1 = gdf_rails.unary_union.intersection(gdf_nhd.unary_union)

        # convert shapely geometries to gpd dataframe
        pts_list = []
        for pt in [pts_0, pts_1]:
            # if geometry is a single point:
            if pt.geom_type == "Point":
                x_coord, y_coord = pt.coords[0]
                p = [
                    {
                        "properties": {"pt_id": 0},
                        "geometry": mapping(Point(x_coord, y_coord)),
                    }
                ]
            else:
                p = [
                    {"properties": {"pt_id": x}, "geometry": mapping(Point(i.x, i.y))}
                    for x, i in enumerate(pt.geoms)
                ]
            p_gdf = gpd.GeoDataFrame.from_features(p)
            pts_list.append(p_gdf)

        # merge both points gdfs, update crs and write out  files
        points = pd.concat(pts_list).pipe(gpd.GeoDataFrame)
        points.crs = gdf_nhd.crs
        points.to_file(str(x_sect_pts))

        gdf_x_pts = points

        # Buffer:
        gdf_nhd["geometry"] = gdf_nhd["geometry"].buffer(25)  # projected units (m)
        gdf_roads["geometry"] = gdf_roads["geometry"].buffer(50)  # projected units (m)
        gdf_rails["geometry"] = gdf_rails["geometry"].buffer(50)  # projected units (m)
        gdf_x_pts["geometry"] = gdf_x_pts["geometry"].buffer(50)  # projected units (m)

        # check to see if streams and roads intersect with valid results
        try:
            # Intersect nhd x roads and nhd x rails
            nhd_x_roads = gpd.overlay(
                gdf_nhd, gdf_roads, how="intersection"
            )  # nhd x roads
            nhd_x_rails = gpd.overlay(
                gdf_nhd, gdf_rails, how="intersection"
            )  # nhd x roads
        except KeyError:
            pass

        try:
            buff_xs_x_roads = gpd.overlay(
                gdf_x_pts, nhd_x_roads, how="intersection"
            )  # nhd x roads
            # nhd x roads x rails
            buff_xs_x_rails = gpd.overlay(gdf_x_pts, nhd_x_rails, how="intersection")
            merge_list = [buff_xs_x_roads, buff_xs_x_rails]
        except (KeyError, NameError, AttributeError):
            merge_list = [buff_xs_x_roads]

        gdf_x_sect_polys = pd.concat(merge_list, sort=True).pipe(gpd.GeoDataFrame)
        gdf_x_sect_polys.crs = gdf_nhd.crs

    else:
        ###---If no rail roads within the HUC---###

        # find unary_unions as points for streams x roads
        pts_0 = gdf_roads.unary_union.intersection(gdf_nhd.unary_union)

        # convert shapely geometries to gpd dataframe
        p = [
            {"properties": {"pt_id": x}, "geometry": mapping(Point(i.x, i.y))}
            for x, i in enumerate(pts_0.geoms)
        ]
        p_gdf = gpd.GeoDataFrame.from_features(p)

        # write out point files
        points = p_gdf.copy()
        points.crs = gdf_nhd.crs
        points.to_file(str(x_sect_pts))

        gdf_x_pts = points.copy()

        # Buffer:
        gdf_nhd["geometry"] = gdf_nhd["geometry"].buffer(25)  # projected units (m)
        gdf_roads["geometry"] = gdf_roads["geometry"].buffer(50)  # projected units (m)
        gdf_x_pts["geometry"] = gdf_x_pts["geometry"].buffer(50)  # projected units (m)

        # Intersect nhd x roads
        nhd_x_roads = gpd.overlay(gdf_nhd, gdf_roads, how="intersection")  # nhd x roads

        # Intersect nhd x roads by buffered intersection points
        buff_xs_x_roads = gpd.overlay(
            gdf_x_pts, nhd_x_roads, how="intersection"
        )  # nhd x roads

        gdf_x_sect_polys = buff_xs_x_roads.copy()
        gdf_x_sect_polys.crs = gdf_nhd.crs

    if not gdf_x_sect_polys.empty:
        # Add common ID field
        # gdf_x_sect_polys['id'] = np.arange(gdf_x_sect_polys.shape[0])
        gdf_x_sect_polys["id"] = 1
        try:
            gdf_x_sect_polys = gdf_x_sect_polys.dissolve(by="id")
        except ValueError:
            # if TopologyException then buffer geometry by 0.01m before dissolving
            # projected units (m)
            gdf_x_sect_polys["geometry"] = gdf_x_sect_polys["geometry"].buffer(0.01)
            gdf_x_sect_polys = gdf_x_sect_polys.dissolve(by="id")

        gdf_x_sect_polys.to_file(x_sect_polys)

        """ passing circular kernel as footprint arg slows down:
                scipy.ndimge.minimum_filter(), so pass 150,150
                square window dims"""
        # Apply the filter:

        # Get the DEM:
        with rasterio.open(str(str_dem_path)) as ds_dem:
            profile = ds_dem.profile.copy()
            arr_dem = ds_dem.read(1)

        # Get the nodata val mask:
        mask = arr_dem == ds_dem.nodata

        # Assign nodata vals to some huge number:
        arr_dem[mask] = 999999.0

        # apply local minima filter
        arr_min = sc.minimum_filter(arr_dem, size=(150, 150))  # footprint=kernel

        # Re-assign nodata vals:
        arr_min[mask] = ds_dem.nodata

        # Write out min. filtered DEM
        with rasterio.open(ds_min_filter, "w", **profile) as dst:
            dst.write_band(1, arr_min)

        # clip ds_min_filter.tif by x-section polys
        WBT.clip_raster_to_polygon(
            ds_min_filter, x_sect_polys, ds_min_clip, maintain_dimensions=True
        )

        # Read DEMs
        dem_min = open_memory_tif(rasterio.open(ds_min_clip, "r").read(1), profile)
        base_dem = open_memory_tif(
            rasterio.open(str(str_dem_path), "r").read(1), profile
        )
        with rasterio.open(dem_merge, "w", **profile) as dst:
            # merge happens in sequence of rasters
            arr_merge, arr_trans = merge.merge([dem_min, base_dem])
            dst.write(arr_merge)

        # logger
        end = round((timer() - st) / 60.0, 2)
        logger.info(
            f"Pre-breach DEM successfully conditioned. Time elapsed: {end} mins"
        )

        return dem_merge


# ===============================================================================
#  Create weight file for TauDEM D8 FAC
# ===============================================================================
def create_wg_from_streamlines(str_streamlines_path, str_dem_path, str_danglepts_path):
    print("Creating weight grid from streamlines:")

    lst_coords = []
    lst_pts = []
    lst_x = []
    lst_y = []

    with fiona.open(str(str_streamlines_path)) as lines:

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
    with rasterio.open(str_dem_path) as ds_dem:
        out_meta = ds_dem.meta.copy()
        out_meta.update(compress="lzw")
        out_meta.update(dtype=rasterio.int16)
        out_meta.update(nodata=-9999)
        out_meta.update(crs=lines.crs)  # shouldn't be necessary

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
    with rasterio.open(str_danglepts_path, "w", **out_meta) as dest:
        dest.write(arr_danglepts, indexes=1)

    return


# ===============================================================================
#  Mega-function for processing a raw DEM
#   1. Breaching and filling
#   2. TauDEM functions
# ===============================================================================
def preprocess_dem(
    root,
    str_streamlines_path,
    dst_crs,
    run_wg,
    run_taudem,
    physio,
    hucID,
    breach_filepath,
    inputProc,
    logger,
):
    try:
        st = timer()
        # << Define all filenames here >>
        str_dem_path = str(root / f"{hucID}_dem_proj.tif")
        str_danglepts_path = str(root / f"{hucID}_wg.tif")
        p = str(root / f"{hucID}_breach_p.tif")
        sd8 = str(root / f"{hucID}_breach_sd8.tif")
        ad8_wg = str(root / f"{hucID}_breach_ad8_wg.tif")
        ad8_no_wg = str(root / f"{hucID}_breach_ad8_no_wg.tif")
        ord_g = str(root / f"{hucID}_breach_ord_g.tif")
        tree = str(root / f"{hucID}_breach_tree")
        coord = str(root / f"{hucID}_breach_coord")
        net = str(root / f"{hucID}_network.shp")
        w = str(root / f"{hucID}_breach_w.tif")
        w_shp = str(root / f"{hucID}_breach_w.shp")
        slp = str(root / f"{hucID}_breach_slp.tif")
        ang = str(root / f"{hucID}_breach_ang.tif")
        dd = str(root / f"{hucID}_hand.tif")
        wshed_physio = str(root / f"{hucID}_breach_w_diss_physio.shp")

        # for debugging
        # print ("01", str_dem_path)
        # print ("02", breach_filepath_tif_tmp)
        # print ("03", breach_filepath_tif_proj)
        # print ("04", breach_filepath_dep)
        # print ("05", str_danglepts_path)
        # print ("11", p)
        # print ("12", sd8)
        # print ("13", ad8_wg)
        # print ("15", ad8_no_wg)
        # print ("16", ord_g)
        # print ("17", tree)
        # print ("18", coord)
        # print ("19", net)
        # print ("20", w)
        # print ("21", slp)
        # print ("22", ang)
        # print ("23", dd)

        """
        This tool is used to remove the sinks (i.e. topographic depressions and flat areas) from
            digital elevation models (DEMs) using a highly efficient and flexible breaching,
            or carving, method.
        Arg Name: InputDEM, type: string, Description: The input DEM name with file extension
        Arg Name: OutputFile, type: string, Description: The output filename with file extension
        Arg Name: MaxDepth, type: float64, Description: The maximum breach channel depth (-1 to ignore)
        Arg Name: MaxLength, type: int, Description: The maximum length of a breach channel (-1 to ignore)
        Arg Name: ConstrainedBreaching, type: bool, Description: Use constrained breaching
        Arg Name: SubsequentFilling, type: bool, Description: Perform post-breach filling
        """

        if run_wg:
            create_wg_from_streamlines(
                str_streamlines_path, str_dem_path, str_danglepts_path
            )

        if run_taudem:
            # ==============  << 2. D8 FDR with TauDEM >> ================       YES
            d8_flow_dir = f'mpiexec -n {inputProc} d8flowdir -fel "{breach_filepath}" -p "{p}" -sd8 "{sd8}"'
            # Submit command to operating system
            logger.info("Running TauDEM D8 Flow Direction...")
            run_tauDEM(d8_flow_dir, logger)
            logger.info("D8 Flow Direction successfully completed")

            # ============= << 3.a AD8 with weight grid >> ================        YES
            # flow accumulation with NHD end nodes is used to derive stream network
            d8_flow_acc_w_grid = f'mpiexec -n {inputProc} AreaD8 -p "{p}" -ad8 "{ad8_wg}" -wg "{str_danglepts_path}" -nc'
            # Submit command to operating system
            logger.info("Running TauDEM D8 FAC (with weight grid)...")
            run_tauDEM(d8_flow_acc_w_grid, logger)
            logger.info("D8 FAC (with weight grid) successfully completed")

            # ============= << 3.b AD8 no weight grid >> ================
            # flow accumulation with-OUT NHD end nodes is used to derive sub-watersheds
            d8_flow_acc_wo_grid = (
                f'mpiexec -n {inputProc} AreaD8 -p "{p}" -ad8 "{ad8_no_wg}" -nc'
            )
            # Submit command to operating system
            logger.info("Running TauDEM D8 FAC (no weights)...")
            run_tauDEM(d8_flow_acc_wo_grid, logger)
            logger.info("D8 FAC (no weights) successfully completed")

            # ============= << 4 StreamReachandWatershed with TauDEM >> ================
            reach_and_watershed = (
                f'mpiexec -n {inputProc} StreamNet -fel "{breach_filepath}" -p "{p}" '
                f'-ad8 "{ad8_no_wg}" -src "{ad8_wg}" -ord "{ord_g}" -tree "{tree}" -coord "{coord}" -net "{net}" -w "{w}"'
            )
            # Submit command to operating system
            logger.info("Running TauDEM Stream Reach and Watershed...")
            run_tauDEM(reach_and_watershed, logger)
            logger.info("Stream Reach and Watershed successfully completed")

            # ============= << 5. Dinf with TauDEM >> =============        YES
            dInf_flow_dir = f'mpiexec -n {inputProc} DinfFlowDir -fel "{breach_filepath}" -ang "{ang}" -slp "{slp}"'
            # Submit command to operating system
            logger.info("Running TauDEM Dinfinity...")
            run_tauDEM(dInf_flow_dir, logger)
            logger.info("Dinfinity successfully completed")

            # ============= << 6. DinfDistanceDown (HAND) with TauDEM >> ============= YES
            # Use original DEM here
            distmeth = "v"
            statmeth = "ave"
            dInf_dist_down = f'mpiexec -n {inputProc} DinfDistDown -fel "{str_dem_path}" -ang "{ang}" -src "{ad8_wg}" -dd "{dd}" -m {statmeth} {distmeth}'
            logger.info("Running TauDEM Dinf Distance Down...")
            # Submit command to operating system
            run_tauDEM(dInf_dist_down, logger)
            logger.info("Dinf Distance Down successfully completed")

            # polygonize watersheds
            watershed_polygonize(w, w_shp, dst_crs, logger)

            # create watershed polys with physiographic attrs
            join_watershed_attrs(w_shp, physio, net, wshed_physio, logger)

            end = round((timer() - st) / 60.0, 2)
            logger.info(
                f"TauDEM preprocessing steps complete. Time elapsed: {end} mins"
            )

    except:
        logger.critical("Unexpected error:", sys.exc_info()[0])
        raise

    return
