import rasterio
import numpy as np
import fiona
from shapely.geometry import LineString, mapping


def gauss_kern(sigma):
    """
    Returns a normalized 2D gauss kernel array for convolutions

     For wavelet curvature calculation
     Chandana Gangodagame
     Wavelet-Compressed Representation of Landscapes for Hydrologic and Geomorphologic Applications
     March 2016IEEE Geoscience and Remote Sensing Letters 13(4):1-6
     DOI: 10.1109/LGRS.2015.2513011

    Args:
        sigma:

    Returns:
    """

    sigma = int(sigma)

    x, y = np.mgrid[-5 * sigma : 5 * sigma, -5 * sigma : 5 * sigma]

    g2x = (
        (1 - x**2 / sigma**2)
        * np.exp(-(x**2 + y**2) / 2 / sigma**2)
        * 1
        / np.sqrt(2 * np.pi * sigma**2)
        / 4
        / sigma
        * np.exp(float(0.5))
    )
    g2y = (
        (1 - y**2 / sigma**2)
        * np.exp(-(x**2 + y**2) / 2 / sigma**2)
        * 1
        / np.sqrt(2 * np.pi * sigma**2)
        / 4
        / sigma
        * np.exp(float(0.5))
    )

    return g2x, g2y


def channel_width_from_bank_pixels(
    df_coords,
    str_streamlines_path,
    str_bankpixels_path,
    str_reachid,
    i_step,
    max_buff,
    str_chanmet_segs,
    logger,
):
    """
    Calculates channel width and sinuosity using parallel offset buffering

    Args:
        df_coords:
        str_streamlines_path:
        str_bankpixels_path:
        str_reachid:
        i_step:
        max_buff:
        str_chanmet_segs:
        logger:

    Returns:
    """

    logger.info("Channel width from bank pixels -- segmented reaches:")

    j = 0
    gp_coords = df_coords.groupby("linkno")

    # Schema for the output properties file:
    schema_output = {
        "geometry": "LineString",
        "properties": {
            "linkno": "int",
            "ch_wid_total": "float",
            "ch_wid_1": "float",
            "ch_wid_2": "float",
            "dist_sl": "float",
            "dist": "float",
            "sinuosity": "float",
            "order": "int",
        },
    }

    # Access the bank pixel layer:
    # open with share=False for multithreading
    with rasterio.open(str(str_bankpixels_path)) as ds_bankpixels:
        # Successive buffer-mask operations to count bank pixels at certain intervals
        lst_buff = range(int(ds_bankpixels.res[0]), max_buff, int(ds_bankpixels.res[0]))
        # Access the streamlines layer:
        with fiona.open(str(str_streamlines_path), "r") as streamlines:
            # Get the crs:
            streamlines_crs = streamlines.crs
            # Open another file to write the output props:
            with fiona.open(
                str(str_chanmet_segs),
                "w",
                "ESRI Shapefile",
                schema_output,
                streamlines_crs,
            ) as output:
                for i_linkno, df_linkno in gp_coords:
                    j += 1
                    i_linkno = int(i_linkno)
                    logger.info("linkno:  {}".format(i_linkno))

                    # << Analysis by reach segments >>
                    # Set up index array to split up df_linkno into segments
                    # (these dictate the reach segment length):
                    # NOTE:  Reach might not be long enough to break up
                    arr_ind = np.arange(
                        i_step, len(df_linkno.index) + 1, i_step
                    )  # NOTE: Change the step for resolution
                    lst_dfsegs = np.split(df_linkno, arr_ind)

                    for i_seg, df_seg in enumerate(
                        lst_dfsegs
                    ):  # looping over each reach segment

                        order = df_seg.order.max()

                        try:
                            order = int(order)
                        except:
                            order = 1

                        arr_x = df_seg.x.values
                        arr_y = df_seg.y.values

                        try:
                            # Create a line segment from endpts in df_seg:
                            ls = LineString(zip(arr_x, arr_y))
                        except:
                            logger.error(
                                "Cannot create a LineString using these points, skipping"
                            )
                            continue

                        try:
                            # Calculate straight line distance:
                            dist_sl = np.sqrt(
                                (arr_x[0] - arr_x[-1]) ** 2
                                + (arr_y[0] - arr_y[-1]) ** 2
                            )
                        except:
                            logger.warning("Error calculated straight line distance")
                            dist_sl = -9999.0

                        dist = ls.length
                        # ratio of sinuous length to straight line length
                        sinuosity = dist / dist_sl

                        lst_tally = []

                        for buff_dist in lst_buff:

                            try:
                                # Watch out for potential geometry errors here:
                                ls_offset_left = ls.parallel_offset(buff_dist, "left")
                                ls_offset_rt = ls.parallel_offset(buff_dist, "right")
                            except:
                                logger.warning("Error performing offset buffer")

                            # Buffer errors can result from complicated line geometry:
                            try:
                                out_left, out_transform = rasterio.mask.mask(
                                    ds_bankpixels, [mapping(ls_offset_left)], crop=True
                                )
                            except:
                                logger.warning("Left offset error")
                                out_left = np.array([0])

                            try:
                                out_rt, out_transform = rasterio.mask.mask(
                                    ds_bankpixels, [mapping(ls_offset_rt)], crop=True
                                )
                            except:
                                logger.warning("Right offset error")
                                out_rt = np.array([0])

                            num_pixels_left = len(out_left[out_left > 0.0])
                            num_pixels_rt = len(out_rt[out_rt > 0.0])

                            # You want the number of pixels gained by each interval:
                            tpl_out = (
                                i_linkno,
                                buff_dist,
                                num_pixels_left,
                                num_pixels_rt,
                            )
                            lst_tally.append(tpl_out)
                            df_tally = pd.DataFrame(
                                lst_tally,
                                columns=[
                                    "linkno",
                                    "buffer",
                                    "interval_left",
                                    "interval_rt",
                                ],
                            )

                        # Calculate weighted average
                        # Only iterate over the top 3 or 2 (n_top) since distance is favored:
                        weighted_avg_left = 0
                        weighted_avg_rt = 0
                        n_top = 2

                        try:
                            for tpl in (
                                df_tally.nlargest(n_top, "interval_left")
                                .iloc[0:2]
                                .itertuples()
                            ):
                                weighted_avg_left += tpl.buffer * (
                                    np.float(tpl.interval_left)
                                    / np.float(
                                        df_tally.nlargest(n_top, "interval_left")
                                        .iloc[0:2]
                                        .sum()
                                        .interval_left
                                    )
                                )
                        except Exception as e:
                            weighted_avg_left = max_buff
                            logger.warning(
                                "Left width set to max. Exception: {} \n".format(e)
                            )

                        try:
                            for tpl in (
                                df_tally.nlargest(n_top, "interval_rt")
                                .iloc[0:2]
                                .itertuples()
                            ):
                                weighted_avg_rt += tpl.buffer * (
                                    np.float(tpl.interval_rt)
                                    / np.float(
                                        df_tally.nlargest(n_top, "interval_rt")
                                        .iloc[0:2]
                                        .sum()
                                        .interval_rt
                                    )
                                )
                        except Exception as e:
                            weighted_avg_rt = max_buff
                            logger.warning(
                                "Right width set to max. Exception: {} \n".format(e)
                            )

                        # Write to the output shapefile here:
                        output.write(
                            {
                                "properties": {
                                    "linkno": i_linkno,
                                    "ch_wid_total": weighted_avg_left + weighted_avg_rt,
                                    "ch_wid_1": weighted_avg_left,
                                    "ch_wid_2": weighted_avg_rt,
                                    "dist_sl": dist_sl,
                                    "dist": dist,
                                    "sinuosity": sinuosity,
                                    "order": int(order),
                                },
                                "geometry": mapping(ls),
                            }
                        )

    return


def bankpixels_from_curvature_window(
    df_coords, str_dem_path, str_bankpixels_path, cell_size, use_wavelet_method, logger
):
    """
    Searchable window with center pixel defined by get_stream_coords_from_features

    Args:
        df_coords:
        str_dem_path:
        str_bankpixels_path:
        cell_size:
        use_wavelet_method:
        logger:

    Returns:
    """
    logger.info("Bank pixels from curvature windows:")

    # Convert df_coords x-y to row-col via DEM affine
    # Loop over center row-col pairs accessing the window
    # Now loop over the linknos to get access grid by window:

    # << PARAMETERS >>
    cell_size = int(cell_size)

    # 3 m:
    w_height = 20  # number of rows
    w_width = 20  # number of columns
    buff = 3  # number of cells
    curve_thresh = 0.30  # good for 3m DEM

    j = 0

    try:
        # select the scale sigma=1
        sigma = 1.0
        g2x1, g2y1 = gauss_kern(sigma)

        with rasterio.open(str_dem_path) as ds_dem:

            # Transform to pixel space
            df_coords["col"], df_coords["row"] = ~ds_dem.transform * (
                df_coords["x"],
                df_coords["y"],
            )
            df_coords[["row", "col"]] = df_coords[["row", "col"]].astype(np.int32)
            df_coords.drop_duplicates(
                ["col", "row"], inplace=True
            )  # rounding to integer
            # total_len = len(df_coords.index)

            out_meta = ds_dem.meta.copy()
            # no need for float32 for bankpixels to save size of output
            out_meta["dtype"] = rasterio.uint8
            out_meta["compress"] = "lzw"

            arr_bankpts = np.zeros(
                [out_meta["height"], out_meta["width"]], dtype=out_meta["dtype"]
            )

            for tpl_row in df_coords.itertuples():

                if tpl_row.order == 5:
                    w_height = 40  # number of rows
                    w_width = 40  # number of columns
                if tpl_row.order >= 6:
                    w_height = 80
                    w_width = 80

                j += 1

                # logger.info('{} | {} -- {}'.format(tpl_row.linkno, j, total_len))

                row_min = np.int(tpl_row.row - np.int(w_height / 2))
                row_max = np.int(tpl_row.row + np.int(w_height / 2))
                col_min = np.int(tpl_row.col - np.int(w_width / 2))
                col_max = np.int(tpl_row.col + np.int(w_width / 2))

                # Now get the DEM specified by this window as a numpy array:
                w = ds_dem.read(1, window=((row_min, row_max), (col_min, col_max)))

                # Then extract the internal part of the window that contains the rotated window
                w[
                    w > 9999999.0
                ] = 0.0  # NoData values may have been corrupted by preprocessing
                w[w < -9999999.0] = 0.0

                # make sure a window of appropriate size was returned from the DEM
                if np.size(w) > 9:

                    if use_wavelet_method:
                        # === Wavelet Curvature from Chandana ===
                        gradfx1 = signal.convolve2d(
                            w, g2x1, boundary="symm", mode="same"
                        )
                        gradfy1 = signal.convolve2d(
                            w, g2y1, boundary="symm", mode="same"
                        )

                        w_curve = gradfx1 + gradfy1

                        # Pick out bankpts:
                        w_curve[w_curve < np.max(w_curve) * curve_thresh] = 0.0

                    else:
                        # Mean Curvature:
                        Zy, Zx = np.gradient(w, cell_size)
                        Zxy, Zxx = np.gradient(Zx, cell_size)
                        Zyy, _ = np.gradient(Zy, cell_size)

                        try:
                            w_curve = (
                                (Zx**2 + 1) * Zyy
                                - 2 * Zx * Zy * Zxy
                                + (Zy**2 + 1) * Zxx
                            )
                            w_curve = -w_curve / (2 * (Zx**2 + Zy**2 + 1) ** (1.5))
                        except:
                            logger.info(
                                "Error calculating Curvature in window:skipping"
                            )
                            continue

                        w_curve[w_curve < np.max(w_curve) * curve_thresh] = 0.0

                    w_curve[w_curve < -99999999.0] = 0.0
                    w_curve[w_curve > 99999999.0] = 0.0

                    w_curve[w_curve > 0.0] = 1.0

                    # Note:  This assumes that the w_curve window is the specified size,
                    # which is not always the case for edge reaches:
                    # arr_bankpts[
                    #   row_min + buff:row_max - buff, col_min + buff:col_max - buff
                    #   ] = w_curve[buff:w_height - buff, buff:w_width - buff]
                    arr_bankpts[
                        row_min + buff : row_max - buff, col_min + buff : col_max - buff
                    ] = w_curve[buff : w_height - buff, buff : w_width - buff]

                    out_meta["nodata"] = 0.0

            logger.info("Writing bank pixels .tif:")
            with rasterio.open(str_bankpixels_path, "w", **out_meta) as dest:
                dest.write(arr_bankpts.astype(rasterio.uint8), indexes=1)

    except Exception as e:
        logger.info(
            "\r\nError in bankpixels_from_curvature_window. Exception: {} \n".format(e)
        )

    return
