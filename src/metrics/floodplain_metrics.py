import rasterio
import rasterio.mask
import geopandas as gpd
from shapely.geometry import mapping, LineString, Point
import fiona
import numpy as np
import pandas as pd
from scipy.ndimage import label

from src import preprocessing


def fp_metrics_chsegs(str_fim_path, str_chwid, str_chanmet_segs, logger):
    """
    Calculate floodplain metrics from 2D cross-sections

    Args:
        str_fim_path: path to the flood inundation raster
        str_chwid: name of the field in the channel segment layer containing pre-calculated channel width values
        str_chanmet_segs: path to the segmented streamline layer containing pre-calculated channel metrics
         (output form channel width from bank pixel method)
        logger: logger instance

    Returns: None
    """
    # Open the stream network segments layer with channel metrics:
    gdf_segs = gpd.read_file(str(str_chanmet_segs))

    lst_fpwid = []
    lst_fprng = []
    lst_geom = []
    lst_min = []
    lst_max = []
    lst_std = []
    lst_rug = []

    # Keep only valid geometries:
    gdf_segs = gdf_segs[gdf_segs.geometry.is_valid]

    # Open the floodplain layer:
    with rasterio.open(str(str_fim_path)) as ds_fim:
        # Loop over each segment:
        for tpl in gdf_segs.itertuples():

            try:
                # if tpl.Index==931: # FOR TESTING
                #      logger.info('pause')

                # Get Xn length based on stream order:
                p_xnlength, p_fitlength = preprocessing.get_xn_length_by_order(tpl.order, True)

                x, y = zip(*mapping(tpl.geometry)["coordinates"])

                # Get the segment midpoint:
                midpt_x = x[int(len(x) / 2)]
                midpt_y = y[int(len(y) / 2)]

                # Build a 1D cross-section from the end points:
                lst_xy = preprocessing.build_xns(y, x, midpt_x, midpt_y, p_xnlength)

                try:
                    # Turn the cross-section into a linestring:
                    fp_ls = LineString([Point(lst_xy[0]), Point(lst_xy[1])])
                except:
                    logger.info("Error converting Xn endpts to LineString")
                    pass

                # Buffer the cross-section to form a 2D rectangle:
                # about half of the line segment straight line distance
                buff_len = tpl.dist_sl / 1.85
                geom_fpls_buff = fp_ls.buffer(buff_len, cap_style=2)
                xn_buff = mapping(geom_fpls_buff)

                # Mask the fp for each feature:
                w_fim, trans_fim = rasterio.mask.mask(ds_fim, [xn_buff], crop=True)
                w_fim = w_fim[0]
                w_fim = w_fim[w_fim != ds_fim.nodata]

                if w_fim[w_fim > 0].size == 0:
                    lst_fpwid.append(-9999.0)
                    lst_fprng.append(-9999.0)
                    lst_geom.append(-9999.0)
                    lst_min.append(-9999.0)
                    lst_max.append(-9999.0)
                    lst_std.append(-9999.0)
                    lst_rug.append(-9999.0)
                    continue

                # OR, Get indices of FIM pixels and use those for the DEM
                # inds_fim=np.where(w_fim!=ds_fim.nodata)

                fp_fim = w_fim[w_fim > 0]  # Assumes is along zero height pixels
                fp_min = fp_fim.min()
                fp_max = fp_fim.max()
                fp_std = fp_fim.std()

                # << Related to mapping the floodplain based on HAND height >>
                # Count the number of pixels in the buffered Xn:
                num_pixels = w_fim.size

                # Calculate area of FP pixels:
                area_pixels = num_pixels * (ds_fim.res[0] ** 2)  # get grid resolution

                # Calculate width by stretching it along the length of the 2D Xn:
                fp_width = area_pixels / (buff_len * 2)
                #    fp_width=0 # For testing purposes

                # Elevation range using HAND heights:
                try:
                    fp_range = fp_max - fp_min
                except:
                    fp_range = 0
                    pass

                # Subtract channel width from fp width:
                fp_width = fp_width - getattr(tpl, str_chwid)
                # If negative, just set it to zero:
                if fp_width < 0.0:
                    fp_width = 0

                # Try calculating roughness (Planar area vs. actual area):
                fp_rug = rugosity(
                    w_fim, ds_fim.res[0], logger
                )  # returns -9999. if error

                lst_min.append(fp_min)
                lst_max.append(fp_max)
                lst_std.append(fp_std)
                lst_fpwid.append(fp_width)
                lst_fprng.append(fp_range)
                lst_rug.append(fp_rug)
                lst_geom.append(tpl.geometry)
            except Exception as e:
                logger.info(f"Error with segment {tpl.Index}: {str(e)}")
                lst_fpwid.append(-9999.0)
                lst_fprng.append(-9999.0)
                lst_geom.append(-9999.0)
                lst_min.append(-9999.0)
                lst_max.append(-9999.0)
                lst_std.append(-9999.0)
                lst_rug.append(-9999.0)
                continue

    # Re-save the channel metrics shapefile with FP metrics added:
    gdf_segs["fp_width_2d"] = lst_fpwid
    gdf_segs["fp_range_2d"] = lst_fprng
    gdf_segs["fp_min_2d"] = lst_min
    gdf_segs["fp_max_2d"] = lst_max
    gdf_segs["fp_std_2d"] = lst_std
    gdf_segs["fp_rug_2d"] = lst_rug
    gdf_segs.to_file(str(str_chanmet_segs))  # [:-4]+'_TEST.shp')


def fim_hand_poly(
    str_hand_path, str_sheds_path, str_reachid, str_fim_path, str_fim_csv, logger
):
    """
    Delineate a FIM from the HAND grid using depth at each polygon (eg, catchment)

    Args:
        str_hand_path:
        str_sheds_path:
        str_reachid:
        str_fim_path:
        str_fim_csv:
        logger:

    Returns:
    """
    # Open the HAND layer:
    with rasterio.open(str(str_hand_path)) as ds_hand:

        out_meta = ds_hand.meta.copy()
        arr_fim = np.empty(
            [out_meta["height"], out_meta["width"]], dtype=out_meta["dtype"]
        )
        arr_fim[:, :] = out_meta["nodata"]

        lst_h = []
        lst_linkno = []
        lst_prov = []

        # Open the catchment polygon layer:
        with fiona.open(str(str_sheds_path), "r") as sheds:
            for shed in sheds:
                # Get the linkno:
                linkno = shed["properties"]["LINKNO"]
                # Get the Province:
                prov = shed["properties"]["PROVINCE"]
                # Get the Drainage Area in km^2:
                da_km2 = shed["properties"]["DSContArea"] / 1000000

                if prov == "COASTAL PLAIN" and da_km2 >= 3 and da_km2 <= 3000:
                    h = 1.65
                elif prov == "PIEDMONT" and da_km2 >= 3 and da_km2 <= 3000:
                    h = (np.log10(da_km2) * 0.471 + 0.523) ** 2
                elif prov == "VALLEY AND RIDGE" and da_km2 >= 3 and da_km2 <= 3000:
                    h = (np.log10(da_km2) * 0.471 + 0.375) ** 2
                elif prov == "APPALACHIAN PLATEAUS" and da_km2 >= 3 and da_km2 <= 3000:
                    h = (np.log10(da_km2) * 0.471 + 0.041) ** 2
                elif prov == "BLUE RIDGE" and da_km2 >= 3 and da_km2 <= 3000:
                    # place holder hand method for blue ridge
                    # h = (np.log10(da_km2) * 0.471 + 0.041) ** 2
                    h = 1.56
                else:
                    lst_h.append(-9999)
                    lst_linkno.append(linkno)
                    lst_prov.append(prov)
                    continue  # skip this catchment

                lst_h.append(h)
                lst_linkno.append(linkno)
                lst_prov.append(prov)

                try:
                    # Mask the bankpts file for each feature:
                    w, out_transform = rasterio.mask.mask(
                        ds_hand, [shed["geometry"]], crop=True
                    )
                    w[(w > h)] = out_meta["nodata"]  # Assign NoData to everywhere else

                    # Now write out the FIM for this shed:
                    w = w[0]
                    shp = np.shape(w)

                    # window bounds in x-y space (west, south, east, north)
                    bounds = rasterio.transform.array_bounds(
                        shp[0], shp[1], out_transform
                    )

                    col_min, row_min = ~ds_hand.transform * (bounds[0], bounds[3])

                    row_min = np.int(row_min)
                    col_min = np.int(col_min)
                    row_max = np.int(row_min + shp[0])
                    col_max = np.int(col_min + shp[1])

                    arr_w = np.empty(
                        [row_max - row_min, col_max - col_min], dtype=out_meta["dtype"]
                    )
                    arr_w[:, :] = arr_fim[row_min:row_max, col_min:col_max]
                    #
                    inds_lt = np.where(arr_fim[row_min:row_max, col_min:col_max] < w)
                    arr_w[inds_lt] = w[inds_lt]

                    # assign the FIM window for this catchment to the total array
                    arr_fim[row_min:row_max, col_min:col_max] = arr_w
                except:
                    logger.info(
                        f"""
                    WARNING: Problem masking HAND grid using catchment Linkno: {linkno}
                    """
                    )

    # Write out the final FIM grid
    out_meta.update(compress="lzw")
    with rasterio.open(
        str_fim_path, "w", tiled=True, blockxsize=512, blockysize=512, **out_meta
    ) as dest:
        dest.write(arr_fim, indexes=1)

    # Write HAND heights to csv
    df_h = pd.DataFrame({str_reachid: lst_linkno, "prov": lst_prov, "h": lst_h})
    df_h.to_csv(str_fim_csv)

    return


def read_fp_xns_shp_and_get_1D_fp_metrics(
    str_xns_path, str_fp_path, str_dem_path, logger
):
    """
    Get FP metrics from 1D cross-section lines

    1. Read Xn file with geopandas, groupby linkno
    2. Linkno extent window like below using rasterio
    3. For each Xn x-y pair interpolate additional points along the length with shapely
    4. Convert to array space
    5. Sample DEM and fp grids as numpy arrays
    6. Calculate metrics

    Args:
        str_xns_path:
        str_fp_path:
        str_dem_path:
        logger:

    Returns:
    """
    # Depth:
    lst_min_d = []
    lst_max_d = []
    lst_rng_d = []
    lst_mean_d = []
    lst_std_d = []
    lst_sum_d = []
    # Elevation:
    lst_min_e = []
    lst_max_e = []
    lst_rng_e = []
    lst_mean_e = []
    lst_std_e = []
    lst_sum_e = []

    lst_id = []
    lst_width = []
    lst_index = []

    # Read xn file:
    logger.info("Reading Xn file:")
    gdf_xns = gpd.read_file(str(str_xns_path))
    # Groupby linkno:
    gp_xns = gdf_xns.groupby("linkno")

    # Access the floodplain and DEM grids:
    with rasterio.open(str(str_dem_path)) as ds_dem:
        with rasterio.open(str(str_fp_path)) as ds_fp:
            # Loop over the linkno groups:
            for linkno, gdf in gp_xns:

                # Loop over the Xns along this linkno:
                for i, tpl in enumerate(gdf.itertuples()):
                    try:
                        # Xn ID and index for saving:
                        lst_id.append(i)
                        lst_index.append(tpl.Index)
                    except Exception as e:
                        print(f"Error with itertuple: {str(e)}")

                    try:
                        # Mask the floodplain grid with each Xn:
                        w_fp, w_trans = rasterio.mask.mask(
                            ds_fp, [mapping(tpl.geometry)], crop=True
                        )
                        w_fp = w_fp[0]
                        w_fp = w_fp[w_fp != ds_fp.nodata]  # ignore nodata vals

                        num_pixels = w_fp.size  # number of fp pixels along the xn
                        tot_width = (
                            num_pixels * ds_fp.res[0]
                        )  # num pixels times cell resolution
                        lst_width.append(tot_width)
                    except:
                        lst_width.append(-9999.0)

                    try:
                        # Relative elevation (depth) metrics:
                        min_depth = w_fp.min()
                        max_depth = w_fp.max()
                        rng_depth = max_depth - min_depth
                        mean_depth = w_fp.mean()
                        std_depth = w_fp.std()
                        sum_depth = w_fp.sum()
                        lst_min_d.append(min_depth)
                        lst_max_d.append(max_depth)
                        lst_rng_d.append(rng_depth)
                        lst_mean_d.append(mean_depth)
                        lst_std_d.append(std_depth)
                        lst_sum_d.append(sum_depth)
                    except:
                        lst_min_d.append(-9999.0)
                        lst_max_d.append(-9999.0)
                        lst_rng_d.append(-9999.0)
                        lst_mean_d.append(-9999.0)
                        lst_std_d.append(-9999.0)
                        lst_sum_d.append(-9999.0)

                    try:
                        # Also mask the DEM to get absolute elevation metrics:
                        w_dem, w_trans = rasterio.mask.mask(
                            ds_dem, [mapping(tpl.geometry)], crop=True
                        )
                        w_dem = w_dem[0]
                        w_dem = w_dem[w_dem != ds_dem.nodata]

                        # Absolute elevation metrics:
                        min_elev = w_dem.min()
                        max_elev = w_dem.max()
                        rng_elev = max_elev - min_elev
                        mean_elev = w_dem.mean()
                        std_elev = w_dem.std()
                        sum_elev = w_dem.sum()
                        lst_min_e.append(min_elev)
                        lst_max_e.append(max_elev)
                        lst_rng_e.append(rng_elev)
                        lst_mean_e.append(mean_elev)
                        lst_std_e.append(std_elev)
                        lst_sum_e.append(sum_elev)
                    except:
                        lst_min_e.append(-9999.0)
                        lst_max_e.append(-9999.0)
                        lst_rng_e.append(-9999.0)
                        lst_mean_e.append(-9999.0)
                        lst_std_e.append(-9999.0)
                        lst_sum_e.append(-9999.0)

            # Initialize fields:
            gdf_xns["xn_id_1dfp"] = -9999.0
            gdf_xns["totwid_1dfp"] = -9999.0
            # Depth:
            gdf_xns["mindep_1dfp"] = -9999.0
            gdf_xns["maxdep_1dfp"] = -9999.0
            gdf_xns["rngdep_1dfp"] = -9999.0
            gdf_xns["meandep_1dfp"] = -9999.0
            gdf_xns["stddep_1dfp"] = -9999.0
            gdf_xns["sumdep_1dfp"] = -9999.0
            # Elevation:
            gdf_xns["minele_1dfp"] = -9999.0
            gdf_xns["maxele_1dfp"] = -9999.0
            gdf_xns["rngele_1dfp"] = -9999.0
            gdf_xns["meanele_1dfp"] = -9999.0
            gdf_xns["stdele_1dfp"] = -9999.0
            gdf_xns["sumele_1dfp"] = -9999.0

            # Add new values:
            gdf_xns.loc[lst_index, "xn_id_1dfp"] = lst_id
            gdf_xns.loc[lst_index, "totwid_1dfp"] = lst_width
            # Depth:
            gdf_xns.loc[lst_index, "mindep_1dfp"] = lst_min_d
            gdf_xns.loc[lst_index, "maxdep_1dfp"] = lst_max_d
            gdf_xns.loc[lst_index, "rngdep_1dfp"] = lst_rng_d
            gdf_xns.loc[lst_index, "meandep_1dfp"] = lst_mean_d
            gdf_xns.loc[lst_index, "stddep_1dfp"] = lst_std_d
            gdf_xns.loc[lst_index, "sumdep_1dfp"] = lst_sum_d
            # Elevation:
            gdf_xns.loc[lst_index, "minele_1dfp"] = lst_min_e
            gdf_xns.loc[lst_index, "maxele_1dfp"] = lst_max_e
            gdf_xns.loc[lst_index, "rngele_1dfp"] = lst_rng_e
            gdf_xns.loc[lst_index, "meanele_1dfp"] = lst_mean_e
            gdf_xns.loc[lst_index, "stdele_1dfp"] = lst_std_e
            gdf_xns.loc[lst_index, "sumele_1dfp"] = lst_sum_e

            # Save it again:
            gdf_xns.to_file(str(str_xns_path))


def hand_analysis_chsegs(
    str_hand_path, str_chanmet_segs, str_src_path, str_fp_path, str_dem_path, logger
):
    """
    Calculation of floodplain metrics by analyzing HAND using 2D cross-sections and
     discerning between channel pixels and floodplain pixels

    Args:
        str_hand_path: Path to the HAND grid .tif
        str_chanmet_segs: Path to the output of the channel_width_from_bank_pixels() func
        str_src_path: Path to the .tif file where stream reaches correspond to linkno values
        str_fp_path: Path to the floodplain grid .tif
        str_dem_path:
        logger: Logger instance for messaging

    Returns: Writes out additional attributes to the file in str_chanmet_segs
    """
    # Open the stream network segments layer with channel metrics:
    gdf_segs = gpd.read_file(str(str_chanmet_segs))

    lst_linkno = []
    # Channel metrics:
    lst_bnk_ht = []
    lst_chn_wid = []
    lst_chn_shp = []
    lst_geom = []  # for testing/creating a new output file
    # FP metrics:
    lst_fpmin = []
    lst_fpmax = []
    lst_fpstd = []
    lst_fprug = []
    lst_fpwid = []
    lst_fprange = []
    lst_fpmin_e = []
    lst_fpmax_e = []
    lst_fpstd_e = []
    lst_fprange_e = []

    # Open the hand layer:
    with rasterio.open(str(str_hand_path)) as ds_hand:

        out_meta = ds_hand.meta.copy()
        arr_chn = np.empty(
            [out_meta["height"], out_meta["width"]], dtype=out_meta["dtype"]
        )
        arr_chn[:, :] = out_meta["nodata"]

        # Access the src layer for excluding channel pixels:
        with rasterio.open(str(str_src_path)) as ds_src:

            res = ds_hand.res[0]

            # Access the floodplain grid and dem
            with rasterio.open(str(str_fp_path)) as ds_fp:
                with rasterio.open(str(str_dem_path)) as ds_dem:

                    # Loop over each segment:
                    for tpl in gdf_segs.itertuples():
                        try:
                            logger.info(f"\t{tpl.Index}")

                            # Get Xn length based on stream order:
                            p_xnlength, p_fitlength = preprocessing.get_xn_length_by_order(
                                tpl.order, False
                            )

                            x, y = zip(*mapping(tpl.geometry)["coordinates"])

                            # Get the segment midpoint:
                            # Can also use this to identify the channel blob
                            midpt_x = x[
                                int(len(x) / 2)
                            ]  # This isn't actually midpt by distance
                            midpt_y = y[int(len(y) / 2)]

                            # Build a 1D cross-section from the end points:
                            lst_xy = preprocessing.build_xns(y, x, midpt_x, midpt_y, p_xnlength)

                            try:
                                # Turn the cross-section into a linestring:
                                fp_ls = LineString([Point(lst_xy[0]), Point(lst_xy[1])])
                            except Exception as e:
                                logger.info(
                                    f"Error converting Xn endpts to LineString: {str(e)}"
                                )
                                pass

                            # Buffer the cross-section to form a 2D rectangle:
                            # about half of the line segment...
                            # ...straight line distance --> AFFECTS LABELLED ARRAY!
                            buff_len = tpl.dist_sl / 2.5
                            geom_fpls_buff = fp_ls.buffer(buff_len, cap_style=2)
                            xn_buff = mapping(geom_fpls_buff)

                            # Mask the hand grid for each feature:
                            w_hand, trans_hand = rasterio.mask.mask(
                                ds_hand, [xn_buff], crop=True
                            )
                            w_hand = w_hand[0]
                            w_hand[w_hand == ds_hand.nodata] = -9999.0

                            # Set up vertical intervals to...
                            # ...slice using 2D cross-section horizontal plane:
                            w_min = w_hand[w_hand > -9999.0].min()
                            w_max = w_hand.max()
                            arr_slices = np.linspace(w_min, w_max, 50)

                            # Also mask the src layer (1 time)...
                            # ...to get the indices of the raster streamline:
                            w_src, trans_src = rasterio.mask.mask(
                                ds_src, [xn_buff], crop=True
                            )
                            w_src = w_src[0]

                            # Also mask the floodplain grid:
                            w_fp, trans_fp = rasterio.mask.mask(
                                ds_fp, [xn_buff], crop=True
                            )
                            w_fp = w_fp[0]

                            # Aaaand, mask the dem:
                            w_dem, trans_dem = rasterio.mask.mask(
                                ds_dem, [xn_buff], crop=True
                            )
                            w_dem = w_dem[0]

                            # Channel pixel indices:
                            src_inds = np.where(w_src == tpl.linkno)
                            # Convert to a set to keep only unique values:
                            src_inds = set(zip(src_inds[0], src_inds[1]))

                            # To produce labeled array:
                            s = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]

                            # Consider the blobs of hand pixels in each slice:
                            lst_count = []
                            lst_width = []
                            # set_inds=set([])
                            lst_inds = []
                            lst_height = []

                            # Begin vertical slicing using 2D cross-section horizontal plane:
                            w_inds = np.indices(w_hand.shape)
                            for i, i_step in enumerate(
                                arr_slices[1:]
                            ):  # skip the first entry

                                # Create a binary array where within step height threshold:
                                w_step = w_hand.copy()
                                w_step[(w_step < i_step) & (w_step > -9999.0)] = 1
                                w_step[w_step != 1] = 0

                                # scipy labeled array:
                                labeled_arr, num_feats = label(w_step, structure=s)

                                # You need to loop over num_feats here and do the test:
                                for feat in np.arange(0, num_feats):
                                    # Get the window indices of each feature:
                                    inds = set(
                                        zip(
                                            w_inds[0][labeled_arr == feat + 1],
                                            w_inds[1][labeled_arr == feat + 1],
                                        )
                                    )

                                    # if they share indices...
                                    # ...consider this blob connected to the channel
                                    if len(src_inds.intersection(inds)) > 0:
                                        lst_count.append(len(inds))
                                        lst_width.append(
                                            len(inds) * (res**2) / tpl.dist_sl
                                        )
                                        lst_height.append(i_step)
                                        lst_inds.append(inds)

                            # End slices here
                            df_steps = pd.DataFrame(
                                {
                                    "count": lst_count,
                                    "height": lst_height,
                                    "width": lst_width,
                                    "inds": lst_inds,
                                }
                            )

                            if len(df_steps.index) < 3:
                                logger.info("Too few slices!")
                                lst_bnk_ht.append(-9999.0)
                                lst_chn_wid.append(-9999.0)
                                lst_chn_shp.append(-9999.0)
                                lst_linkno.append(tpl.linkno)
                                lst_geom.append(tpl.geometry)
                                # FP metrics:
                                lst_fpmax.append(-9999.0)
                                lst_fpmin.append(-9999.0)
                                lst_fpstd.append(-9999.0)
                                lst_fprug.append(-9999.0)
                                lst_fpwid.append(-9999.0)
                                lst_fprange.append(-9999.0)
                                lst_fpmin_e.append(-9999.0)
                                lst_fpmax_e.append(-9999.0)
                                lst_fpstd_e.append(-9999.0)
                                lst_fprange_e.append(-9999.0)
                                continue

                            df_steps["dy"] = df_steps.height.diff()
                            df_steps["dx"] = df_steps.width.diff()
                            df_steps["delta_width"] = df_steps.dx / df_steps.dy
                            indx = df_steps.delta_width.iloc[1:].idxmax() - 1

                            chn_wid = df_steps.width.iloc[indx]
                            bnk_ht = df_steps.height.iloc[indx]
                            chn_shp = np.arctan(
                                bnk_ht / chn_wid
                            )  # sort of entrenchment

                            # Separate the FP and channel pixels using bnk_ht
                            # Channel pixels only:
                            for i_set in df_steps.inds.iloc[0 : indx + 1].tolist():
                                src_inds.update(i_set)

                            # NEED A TUPLE OF 1D ARRAYS FOR ARRAY INDEXING
                            lst1, lst2 = zip(*list(src_inds))

                            # Get the FP pixels without the channel pixels:
                            mask = np.ones_like(w_hand, dtype=bool)
                            mask[lst1, lst2] = False

                            # Relative elevation (HAND):
                            try:
                                w_fp = w_fp[mask]
                                w_fp = w_fp[
                                    w_fp != ds_fp.nodata
                                ]  # also remove nodata vals

                                if w_fp.size == 0:
                                    logger.info("No FP!")
                                    # There's nothing we can do here related to FP:
                                    lst_fpmax.append(-9999.0)
                                    lst_fpmin.append(-9999.0)
                                    lst_fpstd.append(-9999.0)
                                    lst_fprug.append(-9999.0)
                                    lst_fpwid.append(-9999.0)
                                    lst_fprange.append(-9999.0)
                                    lst_fpmin_e.append(-9999.0)
                                    lst_fpmax_e.append(-9999.0)
                                    lst_fpstd_e.append(-9999.0)
                                    lst_fprange_e.append(-9999.0)
                                    # Channel metrics:
                                    lst_bnk_ht.append(bnk_ht)
                                    lst_chn_wid.append(chn_wid)
                                    lst_chn_shp.append(chn_shp)
                                    lst_linkno.append(tpl.linkno)
                                    lst_geom.append(tpl.geometry)
                                    continue

                                # FP width:
                                num_pixels = w_fp.size
                                area_pixels = num_pixels * (
                                    ds_fp.res[0] ** 2
                                )  # get grid resolution
                                # Calculate width by stretching it along the length of the 2D Xn:
                                fp_width = area_pixels / (buff_len * 2)

                                # FP roughness (Planar area vs. actual area):
                                # returns -9999. if error
                                fp_rug = rugosity(w_fp, ds_fp.res[0], logger)

                                # Depth range:
                                fp_max = w_fp.max()
                                fp_min = w_fp.min()
                                fp_std = w_fp.std()
                                fp_range = fp_max - fp_min

                            except Exception as e:
                                logger.error(
                                    f"Error calculating relative elevation FP metrics: {e}"
                                )
                                fp_max = -9999.0
                                fp_min = -9999.0
                                fp_std = -9999.0
                                fp_range = -9999.0
                                fp_rug = -9999.0
                                fp_width = -9999.0

                            try:
                                # Absolute elevation (DEM):
                                w_dem = w_dem[mask]
                                w_dem = w_dem[
                                    w_dem != ds_dem.nodata
                                ]  # also remove nodata vals
                                # Elevation range
                                fp_max_e = w_dem.max()
                                fp_min_e = w_dem.min()
                                fp_std_e = w_dem.std()
                                fp_range_e = fp_max_e - fp_min_e

                            except Exception as e:
                                logger.error(
                                    f"Error calculating absolute elevation FP metrics: {e}"
                                )
                                fp_min_e = -9999.0
                                fp_max_e = -9999.0
                                fp_std_e = -9999.0
                                fp_range_e = -9999.0

                            # Save metrics to lists
                            # Relative elevation:
                            lst_fpmax.append(fp_max)
                            lst_fpmin.append(fp_min)
                            lst_fpstd.append(fp_std)
                            lst_fprug.append(fp_rug)
                            lst_fpwid.append(fp_width)
                            lst_fprange.append(fp_range)
                            # Absolute elevation:
                            lst_fpmin_e.append(fp_min_e)
                            lst_fpmax_e.append(fp_max_e)
                            lst_fpstd_e.append(fp_std_e)
                            lst_fprange_e.append(fp_range_e)
                            # Channel metrics:
                            lst_bnk_ht.append(bnk_ht)
                            lst_chn_wid.append(chn_wid)
                            lst_chn_shp.append(chn_shp)
                            lst_linkno.append(tpl.linkno)
                            lst_geom.append(tpl.geometry)

                            # logger.info('hey')
                        except Exception as e:
                            logger.info(
                                f"Error with segment {tpl.Index}; skipping. {e}"
                            )
                            # sys.exit()
                            # FP metrics:
                            lst_fpmax.append(-9999.0)
                            lst_fpmin.append(-9999.0)
                            lst_fpstd.append(-9999.0)
                            lst_fprug.append(-9999.0)
                            lst_fpwid.append(-9999.0)
                            lst_fprange.append(-9999.0)
                            lst_fpmin_e.append(-9999.0)
                            lst_fpmax_e.append(-9999.0)
                            lst_fpstd_e.append(-9999.0)
                            lst_fprange_e.append(-9999.0)
                            # Channel metrics:
                            lst_bnk_ht.append(-9999.0)
                            lst_chn_wid.append(-9999.0)
                            lst_chn_shp.append(-9999.0)
                            lst_linkno.append(tpl.linkno)
                            lst_geom.append(tpl.geometry)
                            continue

        # Re-save the channel metrics shapefile with FP metrics added:
        gdf_segs["bnk_ht3"] = lst_bnk_ht
        gdf_segs["chn_shp3"] = lst_chn_shp
        gdf_segs["chn_wid3"] = lst_chn_wid
        gdf_segs["fp3_min"] = lst_fpmin
        gdf_segs["fp3_max"] = lst_fpmax
        gdf_segs["fp3_std"] = lst_fpstd
        gdf_segs["fp3_wid"] = lst_fpwid
        gdf_segs["fp3_rug"] = lst_fprug
        gdf_segs["fp3_rng"] = lst_fprange
        gdf_segs["fp3_min_e"] = lst_fpmin_e
        gdf_segs["fp3_max_e"] = lst_fpmax_e
        gdf_segs["fp3_std_e"] = lst_fpstd_e
        gdf_segs["fp3_rng_e"] = lst_fprange_e

        gdf_segs.to_file(str(str_chanmet_segs))


def rugosity(arr, res, logger):
    """
    Actual 3D area divided by 2D planar area gives a measure of terrain complexity or roughness
    Args:
        arr:
        res:
        logger:

    Returns:
    """
    try:
        area3d = (
            (res**2) * (1 + np.gradient(arr) ** 2) ** 0.5
        ).sum()  # actual surface area
        area2d = len(arr) * res**2  # planar surface area
        rug = area3d / area2d
    except:
        logger.info(f"Error in rugosity. arr.shape: {arr.shape}")
        return -9999.0

    return rug
