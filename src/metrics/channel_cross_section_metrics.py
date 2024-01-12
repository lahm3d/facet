from math import ceil, atan

import numpy as np
from osgeo import gdal
import rasterio
import fiona
from fiona.crs import CRS
import pandas as pd
import geopandas as gpd

# from utils import utils


def find_bank_angles(tpl_bfpts, lst_total_slices, xn_len, xn_elev_n, parm_ivert, xn_ptdistance, logger):
    """
    Calculate angle from vertical of left and right banks

    Args:
        tpl_bfpts:
        lst_total_slices:
        xn_len:
        xn_elev_n:
        parm_ivert:
        xn_ptdistance:
        logger:

    Returns:
    """
    try:
        total_slices = len(lst_total_slices)

        # Use second to last slice for bottom estimate
        if total_slices > 2:

            # Interpolate to find left position along bank:
            lf_bottom_ind = lst_total_slices[1][0]

            # LEFT BANK:  Make sure we're within bounds here
            if lf_bottom_ind == 0 or lf_bottom_ind == xn_len:
                lf_angle = 0
            else:
                x1 = lf_bottom_ind - 1
                x2 = lf_bottom_ind
                y1 = xn_elev_n[x1]
                y2 = xn_elev_n[x2]
                yuk = parm_ivert

                lf_bottombank_ind = interp_bank(x1, x2, y1, y2, yuk)

                if abs(lf_bottombank_ind - tpl_bfpts[1]) > 0:
                    # convert radians to degrees
                    lf_angle = (
                        atan(
                            (abs(lf_bottombank_ind - tpl_bfpts[1]))
                            / ((tpl_bfpts[3] - parm_ivert) * xn_ptdistance)
                        )
                        * 57.29578
                    )
                else:
                    lf_angle = 0

            # RIGHT BANK: Interpolate to find left position along bank:
            rt_bottom_ind = lst_total_slices[1][-1]

            # Make sure we're within bounds here
            if rt_bottom_ind == 0 or rt_bottom_ind == xn_len:
                rt_angle = 0
            else:
                x1 = rt_bottom_ind
                x2 = rt_bottom_ind + 1
                y1 = xn_elev_n[x1]
                y2 = xn_elev_n[x2]
                yuk = parm_ivert

                rt_bottombank_ind = interp_bank(x1, x2, y1, y2, yuk)

                if abs(rt_bottombank_ind - tpl_bfpts[2]) > 0:
                    # convert radians to degrees
                    rt_angle = (
                        atan(
                            (abs(rt_bottombank_ind - tpl_bfpts[2]))
                            / ((tpl_bfpts[3] - parm_ivert) * xn_ptdistance)
                        )
                        * 57.29578
                    )
                else:
                    rt_angle = 0

        else:
            # Use bottom slice for bank angle estimate:
            lf_bottom_ind = lst_total_slices[0][0]
            rt_bottom_ind = lst_total_slices[0][-1]

            if abs(lf_bottom_ind - tpl_bfpts[1]) > 0:
                lf_angle = (
                    atan(
                        (abs(lf_bottom_ind - tpl_bfpts[1]) * xn_ptdistance)
                        / tpl_bfpts[3]
                    )
                    * 57.29578
                )  # convert radians to degrees
            else:
                lf_angle = 0

            if abs(rt_bottom_ind - tpl_bfpts[2]) > 0:
                rt_angle = (
                    atan(
                        (abs(rt_bottom_ind - tpl_bfpts[2]) * xn_ptdistance)
                        / tpl_bfpts[3]
                    )
                    * 57.29578
                )  # convert radians to degrees
            else:
                rt_angle = 0

        # NOTE: For now, just set any resulting negative values to -9999,
        # until we figure out what's going on (27Mar2015, SJL)
        if lf_angle < 0:
            lf_angle = -9999.0

        if rt_angle < 0:
            rt_angle = -9999.0

        tpl_angles = (lf_angle, rt_angle)

    except Exception as e:
        logger.info("\r\nError in find_bank_angles. Exception: {} \n".format(e))

    return tpl_angles


def search_right_gt(xnelev, prev_ind, lf_bf):
    """
    Search Xn outward to the right to find the first point greater than the left bank
    Args:
        xnelev:
        prev_ind:
        lf_bf:

    Returns:
    """
    # Search outward from end of previous slice:
    for i in range(prev_ind + 1, len(xnelev), 1):

        if xnelev[i] > lf_bf:
            bank_ind = i
            break
        else:  # flag it so you can skip this one
            bank_ind = -9999

    return bank_ind


def search_left_gt(xnelev, prev_ind, rt_bf):
    """
    Search Xn outward to the left to find the first point greater than the right bank

    Args:
        xnelev:
        prev_ind:
        rt_bf:

    Returns:
    """
    # Search outward from end of previous slice:
    for i in range(prev_ind, 0, -1):

        if xnelev[i] > rt_bf:
            bank_ind = i
            break
        else:  # flag it so you can skip this one
            bank_ind = -9999

    return bank_ind


def interp_bank(x1, x2, y1, y2, y_uk):
    """
    Interpolate to find positions of right/left bank

    Args:
        x1:
        x2:
        y1:
        y2:
        y_uk:

    Returns:
    """
    x_uk = (((x2 - x1) * (y_uk - y1)) / (y2 - y1)) + x1
    return x_uk


def find_bank_ratio_method(lst_total, ratio_threshold, xnelev_zero, slp_thresh, logger):
    """
    Search for banks via a slope threshold and vertical slices. Compares the length of the last gtzero slice
     (num of indices) vs. the previous slice

    Args:
        lst_total: a list of 1D array slice index values
        ratio_threshold:
        xnelev_zero: the Xn elevation profile normalized to zero
        slp_thresh:
        logger:

    Returns: a tuple of bankfull points (left_index, right_index, height)

    """
    tpl_bfpts = ()  # output tuple
    num_slices = (
        len(lst_total) - 1
    )  # total number of slices, each at a height of param_vertstep
    xn_len = len(xnelev_zero) - 1  # length of Xn

    try:
        if num_slices > 2 and len(lst_total[num_slices - 1]) > 2:
            top_area = len(lst_total[num_slices])
            below_area = len(lst_total[num_slices - 1])

            # Check the ratio:
            this_ratio = float(top_area) / float(below_area)

            # USE THIS TO DRIVE THE BANK BREAK DETERMINATION INSTEAD
            if this_ratio > ratio_threshold:
                # Find end indices of this and of previous slice:
                prev_lf_ind = lst_total[num_slices - 1][0]
                prev_rt_ind = lst_total[num_slices - 1][-1]

                this_lf_ind = lst_total[num_slices][0]
                this_rt_ind = lst_total[num_slices][-1]

                # Bottom left and right for searching:
                bottom_lf = lst_total[0][0]
                bottom_rt = lst_total[0][-1]

                # First derivative is slope:
                lf_arr = np.array(xnelev_zero[this_lf_ind : prev_lf_ind + 1])
                rt_arr = np.array(xnelev_zero[prev_rt_ind - 1 : this_rt_ind])
                firstdiff_left = np.diff(lf_arr)
                firstdiff_right = np.diff(rt_arr)

                # Set both indices to negative 1 initially:
                rt_bank_ind = -1
                lf_bank_ind = -1

                # Look for the first occurrence of a very small slope value in both directions:
                for r, this_rt in enumerate(firstdiff_right):
                    if this_rt < slp_thresh:
                        rt_bank_ind = r + prev_rt_ind - 1
                        break

                firstdiff_left_rev = firstdiff_left[::-1]

                for r, this_lf in enumerate(firstdiff_left_rev):
                    if this_lf > -slp_thresh:
                        lf_bank_ind = prev_lf_ind - r
                        break

                # Make sure rt_bank_ind is not greater than total xn length
                if prev_lf_ind > 0 and prev_rt_ind < xn_len:

                    # Find the smallest height of the two:
                    if (
                            rt_bank_ind > 0 > lf_bank_ind
                    ):  # only the right index exists

                        # Interpolate to find left bankfull:
                        bf_height = xnelev_zero[rt_bank_ind]
                        lf_x2 = search_left_gt(xnelev_zero, bottom_rt, bf_height)

                        if lf_x2 != -9999:
                            lf_x1 = lf_x2 + 1

                            lf_y1 = xnelev_zero[lf_x1]
                            lf_y2 = xnelev_zero[lf_x2]

                            if lf_y1 == lf_y2:
                                lfbf_ind = lf_x1
                            else:
                                lfbf_ind = interp_bank(
                                    lf_x1, lf_x2, lf_y1, lf_y2, bf_height
                                )

                            tpl_bfpts = (lfbf_ind, rt_bank_ind, bf_height)

                    elif (
                            lf_bank_ind > 0 > rt_bank_ind
                    ):  # only the left index exists

                        # Interpolate to find right bank index:
                        bf_height = xnelev_zero[lf_bank_ind]
                        rt_x2 = search_right_gt(xnelev_zero, bottom_lf, bf_height)

                        if rt_x2 != -9999:
                            rt_x1 = rt_x2 - 1

                            rt_y1 = xnelev_zero[rt_x1]
                            rt_y2 = xnelev_zero[rt_x2]

                            if rt_y1 == rt_y2:
                                rtbf_ind = rt_x1
                            else:
                                rtbf_ind = interp_bank(
                                    rt_x1, rt_x2, rt_y1, rt_y2, bf_height
                                )

                            tpl_bfpts = (lf_bank_ind, rtbf_ind, bf_height)

                    elif (
                        rt_bank_ind > 0
                        and lf_bank_ind > 0
                        and xnelev_zero[rt_bank_ind] < xnelev_zero[lf_bank_ind]
                    ):  # right is smaller than left

                        # Interpolate to find left bankfull:
                        bf_height = xnelev_zero[rt_bank_ind]
                        # search all the way across
                        lf_x2 = search_left_gt(xnelev_zero, bottom_rt, bf_height)

                        # find the index that's just smaller than bank height on left side FASTER
                        # lf_x2 = search_left_lt(xnelev_zero, lf_bank_ind, bf_height)

                        if lf_x2 != -9999:
                            lf_x1 = lf_x2 + 1

                            lf_y1 = xnelev_zero[lf_x1]
                            lf_y2 = xnelev_zero[lf_x2]

                            if lf_y1 == lf_y2:
                                lfbf_ind = lf_x1
                            else:
                                lfbf_ind = interp_bank(
                                    lf_x1, lf_x2, lf_y1, lf_y2, bf_height
                                )

                            tpl_bfpts = (lfbf_ind, rt_bank_ind, bf_height)

                    elif (
                        rt_bank_ind > 0
                        and lf_bank_ind > 0
                        and xnelev_zero[lf_bank_ind] < xnelev_zero[rt_bank_ind]
                    ):  # left is smaller than right

                        # Interpolate to find right bank index:
                        bf_height = xnelev_zero[lf_bank_ind]
                        rt_x2 = search_right_gt(
                            xnelev_zero, bottom_lf, bf_height
                        )  # Searches all the way across channel

                        if rt_x2 != -9999:
                            rt_x1 = rt_x2 - 1

                            rt_y1 = xnelev_zero[rt_x1]
                            rt_y2 = xnelev_zero[rt_x2]

                            if rt_y1 == rt_y2:
                                rtbf_ind = rt_x1
                            else:
                                rtbf_ind = interp_bank(
                                    rt_x1, rt_x2, rt_y1, rt_y2, bf_height
                                )

                            tpl_bfpts = (lf_bank_ind, rtbf_ind, bf_height)

                    elif (
                        rt_bank_ind > 0
                        and lf_bank_ind > 0
                        and xnelev_zero[lf_bank_ind] == xnelev_zero[rt_bank_ind]
                    ):  # they're exactly equal
                        # logger.info 'they are the same!'
                        bf_height = xnelev_zero[lf_bank_ind]
                        tpl_bfpts = (lf_bank_ind, rt_bank_ind, bf_height)

    except Exception as e:
        logger.info("\r\nError in find_bank_ratio_method. Exception: {} \n".format(e))

    return tpl_bfpts


def is_contiguous(gtzero_inds):
    """
    Check for continuity in vertical cross-section slices. Used by analyze_elev function

    Args:
        gtzero_inds:

    Returns:
    """
    if (np.max(gtzero_inds) - np.min(gtzero_inds)) == np.count_nonzero(gtzero_inds) - 1:
        # Contiguous, continue
        bool_cont = True

    else:
        # Not contiguous, trim off the extras
        bool_cont = False

    return bool_cont


def analyze_xnelev(
    df_xn_elev,
    param_ivert,
    xn_ptdist,
    param_ratiothreshold,
    param_slpthreshold,
    nodata_val,
    logger,
):
    """
    Analyze the elevation profile of each Xn and determine metrics.
         1. Loop over Xn list
         2. Normalize Xn to zero
         3. Loop over vertical slices using a step set by param_ivert
         4. Find contiguous sets of indices for each slice
         5. Search for slope break
         6. If slope break exists, determine "bankfull" locations along Xn

    Args:
        df_xn_elev:
        param_ivert:
        xn_ptdist:
        param_ratiothreshold:
        param_slpthreshold:
        nodata_val:
        logger:

    Returns: List of tuples containing metrics at each cross-section

    tpl_bankfullpts
    lst_total_cnt
    this_linkno
    tpl_bankfullpts[3] + np.min(tpl_row.elev)
    tpl_bankangles
    bf_area
    ch_width
    overbank_ratio
    total_arearatio

    """
    lst_bfmetrics = []  # list of tuples to contain output

    try:
        for tpl_row in df_xn_elev.itertuples():

            this_linkno = tpl_row.linkno
            # this_order = tpl_row.strmord

            # A list to store the total number of indices/blocks in a Xn:
            lst_total_cnt = []

            arr_elev = tpl_row.elev
            arr_elev = arr_elev[arr_elev != np.float32(nodata_val)]

            # Normalize elevation to zero:
            thisxn_norm = tpl_row.elev - np.min(tpl_row.elev)
            # Below is if you're using the breached DEM:
            # if this_order < 5:
            #     thisxn_norm = arr_elev - np.min(arr_elev)
            #     thisxn_norm = tpl_row.elev - np.min(tpl_row.elev)
            # else:
            #     # for order>5:
            #     # (THIS ASSUMES YOU'RE USING THE BREACHED DEM, may not be necessary otherwise)
            #     thisxn_norm = tpl_row.elev - np.partition(tpl_row.elev, 2)[2]
            #     # then any negatives make zero:
            #     thisxn_norm[thisxn_norm<0]=0

            # Loop from zero to max(this_xn_norm) using a pre-defined vertical step (0.2 m):
            for this_slice in np.arange(0.0, np.max(thisxn_norm), param_ivert):

                # The indices of positives:
                gtzero_indices = np.nonzero((this_slice - thisxn_norm) > 0)[
                    0
                ]  # Zero index get the first element of the returned tuple

                # Use funtion to check if contiguous:
                if np.size(gtzero_indices) == 0:  # the first loop only

                    # get the index of the zero value:
                    lst_total_cnt.append(np.where(thisxn_norm == 0)[0])
                    prev_val = lst_total_cnt[0][0]

                elif is_contiguous(gtzero_indices):

                    # Yes, it is contiguous
                    # Save it to the total count:
                    lst_total_cnt.append(gtzero_indices)
                    prev_val = gtzero_indices[
                        0
                    ]  # just need any value from the contiguous array

                else:
                    # No, it's not contiguous
                    # Find the contiguous part of the slice:
                    tpl_parts = np.array_split(
                        gtzero_indices, np.where(np.diff(gtzero_indices) != 1)[0] + 1
                    )  # splits the contiguous elements into separate tuple elements

                    # Find the one that contains an element of the previous slice:
                    # if prev_val in [this_arr for this_arr in tpl_parts]:
                    for this_arr in tpl_parts:
                        if prev_val in this_arr[:]:
                            lst_total_cnt.append(this_arr)
                            prev_val = this_arr[0]
                            break

                tpl_bankfullpts = find_bank_ratio_method(
                    lst_total_cnt,
                    param_ratiothreshold,
                    thisxn_norm,
                    param_slpthreshold,
                    logger,
                )

                if tpl_bankfullpts:
                    # BF points are so close that rounding to int gives identical values
                    if tpl_bankfullpts[0] - tpl_bankfullpts[1] == 0:
                        break

                    # Add Xn number to the output:
                    # normalized elevation profile
                    # xn_elev_norm = tpl_thisxn[2] - np.min(tpl_thisxn[2])
                    xn_length = len(thisxn_norm)

                    # Bank points tuple:
                    #  tpl_bankfullpts = (tpl_thisxn[3],) + tpl_bankfullpts # add local ID
                    tpl_bankfullpts = (tpl_row.Index,) + tpl_bankfullpts

                    # Find bank angles:
                    tpl_bankangles = find_bank_angles(
                        tpl_bankfullpts,
                        lst_total_cnt,
                        xn_length,
                        thisxn_norm,
                        param_ivert,
                        xn_ptdist,
                        logger,
                    )

                    # Estimate bankfull area:
                    # (Bank height - xn_elev_norm[i])*xn_ptdist
                    # Round up on left, round down on right
                    bf_area = 0
                    lst_bf_rng = range(
                        int(ceil(tpl_bankfullpts[1])), int(tpl_bankfullpts[2]) + 1, 1
                    )

                    for i in lst_bf_rng:
                        bf_area += (tpl_bankfullpts[3] - thisxn_norm[i]) * xn_ptdist

                    # Channel width:
                    ch_width = (tpl_bankfullpts[2] - tpl_bankfullpts[1]) * xn_ptdist

                    # Overbank ratio:
                    try:
                        overbank_ratio = len(lst_total_cnt[-1]) / (
                            tpl_bankfullpts[2] - tpl_bankfullpts[1]
                        )
                    except:
                        overbank_ratio = -9999.0

                    if bf_area == 0:
                        bf_area = -9999.0
                        total_arearatio = -9999.0
                    else:
                        # Also try area under entire Xn length relative to BF area:
                        total_xn_area = sum(thisxn_norm * xn_ptdist)
                        try:
                            total_arearatio = (total_xn_area - bf_area) / bf_area
                        except:
                            total_arearatio = -9999.0

                    tpl_metrics = (
                        tpl_bankfullpts
                        + (lst_total_cnt,)
                        + (this_linkno,)
                        + (tpl_bankfullpts[3] + np.min(tpl_row.elev),)
                        + tpl_bankangles
                        + (bf_area,)
                        + (ch_width,)
                        + (overbank_ratio,)
                        + (total_arearatio,)
                    )

                    lst_bfmetrics.append(tpl_metrics)

                    break  # no need to keep slicing here, unless we want to try for FP analysis

    except Exception as e:
        logger.info("\r\nError in analyze_xn_elev. Exception: {} \n".format(e))
        pass

    return lst_bfmetrics


def interpolate(arr_in, ind_val):
    """

    Args:
        arr_in:
        ind_val:

    Returns:
    """
    if ind_val == np.ceil(ind_val):
        out_val = arr_in[int(np.ceil(ind_val))]
    else:
        # it will always be divided by 1
        out_val = arr_in[int(ind_val)] + (ind_val - int(ind_val)) * (
            arr_in[int(np.ceil(ind_val))] - arr_in[int(ind_val)]
        )

    return out_val


def chanmetrics_bankpts(
    df_xn_elev,
    channel_xns,
    dem,
    bank_points,
    parm_ivert,
    XnPtDist,
    parm_ratiothresh,
    parm_slpthresh,
    epsg,
    logger,
):
    """
    Calculate channel metrics based on the bankpoint slope-threshold method at each Xn,
     writing the bank points to a shapefile

    Args:
        df_xn_elev:
        channel_xns:
        dem:
        bank_points:
        parm_ivert:
        XnPtDist:
        parm_ratiothresh:
        parm_slpthresh:
        logger:
        epsg:

    Returns:
    """

    logger.info("Channel metrics from bank points:")

    # << BEGIN LOOP >>
    # Do the rest by looping in strides, rather than all at once, to conserve memory:
    # (possibly using multiprocessing)
    xn_count = gpd.read_file(channel_xns).shape[0]

    # Striding:
    arr_strides = np.linspace(0, xn_count, int(xn_count / 100))
    arr_strides = np.delete(arr_strides, 0)

    # Now loop over the linknos to get access grid by window:
    with rasterio.open(dem) as ds_dem:

        nodata_val = ds_dem.nodata

        # Define the schema for the output bank points shapefile:
        properties_dtypes = {
            "xn_num": "int",
            "linkno": "int",
            "bank_hght": "float",
            "bank_elev": "float",
            "bnk_ang_1": "float",
            "bnk_ang_2": "float",
            "bf_area": "float",
            "chan_width": "float",
            "obank_rat": "float",
            "area_ratio": "float",
        }
        schema = {"geometry": "Point", "properties": properties_dtypes}

        with fiona.open(
            bank_points,
            "w",
            driver="ESRI Shapefile",
            crs=CRS.from_epsg(int(epsg)),
            schema=schema,
        ) as bankpts:
            j = 0
            for indx in arr_strides:
                df_xn_elev_n = df_xn_elev.iloc[j : int(indx)]
                j = int(indx) + 1
                logger.info("\tIndex {} - {}/{}".format(j, int(indx), xn_count))

                # << INTERPOLATE XNs >>
                interpolate_columns = [
                    "xn_no",
                    "left_ind",
                    "right_ind",
                    "bank_height",
                    "slices",
                    "linkno",
                    "bank_elev",
                    "lf_bank_ang",
                    "rt_bank_ang",
                    "bankful_area",
                    "chan_width",
                    "overbank_ratio",
                    "area_ratio",
                ]

                df_bank_metrics = pd.DataFrame(
                    analyze_xnelev(
                        df_xn_elev_n,
                        parm_ivert,
                        XnPtDist,
                        parm_ratiothresh,
                        parm_slpthresh,
                        nodata_val,
                        logger,
                    ),
                    columns=interpolate_columns,
                )

                df_bank_metrics.set_index("xn_no", inplace=True)

                df_map = pd.merge(
                    df_xn_elev, df_bank_metrics, left_index=True, right_index=True
                )

                lst_lfbank_row = []
                lst_lfbank_col = []
                lst_rtbank_row = []
                lst_rtbank_col = []

                for tpl_row in df_map.itertuples():
                    lst_lfbank_row.append(interpolate(tpl_row.xn_row, tpl_row.left_ind))
                    lst_lfbank_col.append(interpolate(tpl_row.xn_col, tpl_row.left_ind))
                    lst_rtbank_row.append(
                        interpolate(tpl_row.xn_row, tpl_row.right_ind)
                    )
                    lst_rtbank_col.append(
                        interpolate(tpl_row.xn_col, tpl_row.right_ind)
                    )

                df_map["lfbank_row"] = pd.Series(lst_lfbank_row).values
                df_map["lfbank_col"] = pd.Series(lst_lfbank_col).values
                df_map["rtbank_row"] = pd.Series(lst_rtbank_row).values
                df_map["rtbank_col"] = pd.Series(lst_rtbank_col).values

                # Transform to pixel space:
                df_map["lfbank_x"], df_map["lfbank_y"] = ds_dem.transform * (
                    df_map["lfbank_col"],
                    df_map["lfbank_row"],
                )
                df_map["rtbank_x"], df_map["rtbank_y"] = ds_dem.transform * (
                    df_map["rtbank_col"],
                    df_map["rtbank_row"],
                )

                for tpl_row in df_map.itertuples():

                    tpl_left = (tpl_row.lfbank_x, tpl_row.lfbank_y)
                    tpl_right = (tpl_row.rtbank_x, tpl_row.rtbank_y)

                    # the shapefile geometry use (lon,lat) Requires a list of x-y tuples
                    lf_pt = {"type": "Point", "coordinates": tpl_left}
                    rt_pt = {"type": "Point", "coordinates": tpl_right}

                    prop_lf = {
                        "xn_num": int(tpl_row.Index),
                        "linkno": int(tpl_row.linkno_x),
                        "bank_hght": tpl_row.bank_height,
                        "bank_elev": tpl_row.bank_elev,
                        "bnk_ang_1": tpl_row.lf_bank_ang,
                        "bf_area": tpl_row.bankful_area,
                        "bnk_ang_2": -9999.0,
                        "chan_width": tpl_row.chan_width,
                        "obank_rat": tpl_row.overbank_ratio,
                        "area_ratio": tpl_row.area_ratio,
                    }

                    prop_rt = {
                        "xn_num": int(tpl_row.Index),
                        "linkno": int(tpl_row.linkno_x),
                        "bank_hght": tpl_row.bank_height,
                        "bank_elev": tpl_row.bank_elev,
                        "bnk_ang_2": tpl_row.rt_bank_ang,
                        "bf_area": tpl_row.bankful_area,
                        "bnk_ang_1": -9999.0,
                        "chan_width": tpl_row.chan_width,
                        "obank_rat": tpl_row.overbank_ratio,
                        "area_ratio": tpl_row.area_ratio,
                    }

                    bankpts.write({"geometry": lf_pt, "properties": prop_lf})
                    bankpts.write({"geometry": rt_pt, "properties": prop_rt})


def read_xns_shp_and_get_dem_window(channel_xns, dem, logger):
    """
    Read an existing Xn file, calculate xy bounds for each linkno and read the DEM
     according to that window

    Args:
        channel_xns:
        dem:
        logger:

    Returns:
    """
    min_nodata_thresh = -99999.0
    max_nodata_thresh = 99999.0

    logger.info("Reading and interpolating elevation along Xn's:")

    lst_linknos = []
    lst_x1 = []
    lst_y1 = []
    lst_x2 = []
    lst_y2 = []
    lst_strmord = []

    #    start_time = timeit.default_timer()
    # First get all linknos:
    with fiona.open(channel_xns, "r") as xn_shp:
        # Read each feature line:
        for line in xn_shp:
            lst_linknos.append(line["properties"]["linkno"])
            lst_x1.append(line["geometry"]["coordinates"][0][0])
            lst_y1.append(line["geometry"]["coordinates"][0][1])
            lst_x2.append(line["geometry"]["coordinates"][1][0])
            lst_y2.append(line["geometry"]["coordinates"][1][1])
            lst_strmord.append(line["properties"]["strmord"])

    df_coords = pd.DataFrame(
        {
            "linkno": lst_linknos,
            "x1": lst_x1,
            "y1": lst_y1,
            "x2": lst_x2,
            "y2": lst_y2,
            "strmord": lst_strmord,
        }
    )

    # Now loop over the linknos to get access grid by window:
    with rasterio.open(dem) as ds_dem:

        nodata_val = (
            ds_dem.nodata
        )  # NODATA val must be defined for this to return anything

        # Get bounds of DEM (left, bottom, right, top):
        bnds = ds_dem.bounds

        # Check the min and max of the coordinates in df_coords-
        # -and remove any cross-sections that extend beyond DEM:
        df_coords["min_x"] = df_coords[["x1", "x2"]].min(axis=1)
        df_coords["max_x"] = df_coords[["x1", "x2"]].max(axis=1)
        df_coords["min_y"] = df_coords[["y1", "y2"]].min(axis=1)
        df_coords["max_y"] = df_coords[["y1", "y2"]].max(axis=1)

        # check min/max_x against bnds[0] and bnds[2] and min/max_y against bnds[1] and bnds[3]
        # min_x > bnds[0], max_x < bnds[2], min_y > bnds[1], max_y < bnds[3]

        df_coords = df_coords[
            (df_coords["min_x"] > bnds[0])
            & (df_coords["max_x"] < bnds[2])
            & (df_coords["min_y"] > bnds[1])
            & (df_coords["max_y"] < bnds[3])
        ]
        # clean columns
        df_coords = df_coords.drop(["min_x", "max_x", "min_y", "max_y"], axis=1)

        # Transform to pixel space
        df_coords["col1"], df_coords["row1"] = ~ds_dem.transform * (
            df_coords["x1"],
            df_coords["y1"],
        )
        df_coords["col2"], df_coords["row2"] = ~ds_dem.transform * (
            df_coords["x2"],
            df_coords["y2"],
        )

        ## OR:
        gp_coords = df_coords.groupby("linkno")

        lst_all_zi = []
        j = 0

        for linkno, df_linkno in gp_coords:
            row_min = int(df_linkno[["row1", "row2"]].min(axis=0).min())
            row_max = int(df_linkno[["row1", "row2"]].max(axis=0).max())
            col_min = int(df_linkno[["col1", "col2"]].min(axis=0).min())
            col_max = int(df_linkno[["col1", "col2"]].max(axis=0).max())
            strmord = int(df_linkno.strmord.iloc[0])

            # Now get the DEM specified by this window as a numpy array:
            w = ds_dem.read(1, window=((row_min, row_max + 1), (col_min, col_max + 1)))

            w_min = np.min(w)
            w_max = np.max(w)

            if w_min < min_nodata_thresh:
                nodata_val = w_min
            elif w_max > max_nodata_thresh:
                nodata_val = w_max

            # NOW loop over each Xn:
            for tpl_xn in df_linkno.itertuples():
                j += 1
                xn_len = int(
                    np.hypot(tpl_xn.col2 - tpl_xn.col1, tpl_xn.row2 - tpl_xn.row1)
                )
                lst_xnrow = np.linspace(
                    tpl_xn.row1 - row_min, tpl_xn.row2 - row_min, xn_len
                )
                lst_xncol = np.linspace(
                    tpl_xn.col1 - col_min, tpl_xn.col2 - col_min, xn_len
                )

                # this is always 1 cell or equivalent to cell_size in meters/feet
                # xnptdist = xn_len/len(lst_xnrow)
                try:
                    arr_zi = w[
                        lst_xnrow.astype(int), lst_xncol.astype(int)
                    ]  # nearest-neighbor
                except:
                    continue

                # Remove possible no data values:NOTE:  They may not be defined in the original file
                arr_zi = arr_zi[arr_zi != np.float32(nodata_val)]

                # if it only has less than 5 elevation measurements along this Xn, skip it
                if arr_zi.size < 5:
                    continue

                # Convert these from window row/col to raster row/col for bankpt use:
                for i, xnrow in enumerate(lst_xnrow):
                    lst_xnrow[i] = lst_xnrow[i] + row_min
                    lst_xncol[i] = lst_xncol[i] + col_min

                tpl_out = (linkno, arr_zi, lst_xnrow, lst_xncol, strmord)
                lst_all_zi.append(tpl_out)

    # print('\tTotal Xn\'s:  {}'.format(i))
    # print('\tTime interpolating elevation along Xn\'s:'+ str(timeit.default_timer()-start_time))

    return pd.DataFrame(
        lst_all_zi, columns=["linkno", "elev", "xn_row", "xn_col", "strmord"]
    )


def derive(channel_xns, dem, bank_points, params, epsg, logger):
    XnPtDist = 1 # same as resolution

    # channel_xns, dem, bank_points, params, epsg = Paths.channel_xns, Paths.dem, Paths.bank_points, Config.methods['cross_section'], Config.spatial_ref['epsg']

    df_xn_elev = read_xns_shp_and_get_dem_window(channel_xns, dem, logger)
    # df_xn_elev.to_csv(str_csv_path)


    chanmetrics_bankpts(
        df_xn_elev,
        channel_xns,
        dem,
        bank_points,
        params['parm_ivert'],
        XnPtDist,
        params['parm_ratiothresh'],
        params['parm_slpthresh'],
        epsg,
        logger,
        )
