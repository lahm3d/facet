import sys
from math import isinf, sqrt
from timeit import default_timer as timer

import fiona
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString



def get_stream_coords_from_features(network, cell_size, reach_id, order_id, csv_output):
    """


    Args:
        network:
        cell_size:
        reach_id:
        order_id:

    Returns:
    """

    if csv_output.is_file():
        return pd.read_csv(csv_output)
    else:
        lst_df_final = []

        p_interp_spacing = int(
            cell_size
        )  # 3 # larger numbers would simulate a more smoothed reach
        j = 0  # prog bar

        # Open the streamlines shapefile:
        with fiona.open(network, "r") as streamlines:
            tot = len(streamlines)
            for line in streamlines:
                j += 1
                line_shply = LineString(line["geometry"]["coordinates"])

                length = line_shply.length  # units depend on crs

                if length > 9:  # Skip small ones. NOTE: This value is dependent on CRS!!

                    i_linkno = line["properties"][reach_id]
                    i_order = line["properties"][order_id]

                    # Smoothing reaches via Shapely:
                    if i_order <= 3:
                        line_shply = line_shply.simplify(5.0, preserve_topology=False)
                    elif i_order == 4:
                        line_shply = line_shply.simplify(10.0, preserve_topology=False)
                    elif i_order == 5:
                        line_shply = line_shply.simplify(20.0, preserve_topology=False)
                    elif i_order >= 6:
                        line_shply = line_shply.simplify(30.0, preserve_topology=False)

                    length = line_shply.length

                    # p_interp_spacing in projection units
                    int_pts = np.arange(0, length, p_interp_spacing)
                    # print(i_linkno, int_pts.shape)

                    lst_x = []
                    lst_y = []
                    lst_linkno = []
                    lst_order = []
                    for i in int_pts:
                        i_pt = np.array(line_shply.interpolate(i)).item()
                        lst_x.append(i_pt.x)
                        lst_y.append(i_pt.y)
                        lst_linkno.append(i_linkno)
                        lst_order.append(i_order)

                    df_coords = pd.DataFrame(
                        {"x": lst_x, "y": lst_y, "linkno": lst_linkno, "order": lst_order}
                    )
                    # potential duplicates due to interpolation
                    df_coords.drop_duplicates(subset=["x", "y"], inplace=True)
                    lst_df_final.append(df_coords)

            df_final = pd.concat(lst_df_final)
            df_final.to_csv(csv_output, index=False)

        return df_final


def build_xns(lstThisSegmentRows, lstThisSegmentCols, midPtCol, midPtRow, p_xnlength):
    """
    Build cross-sections based on vector features

    Args:
        lstThisSegmentRows:
        lstThisSegmentCols:
        midPtCol:
        midPtRow:
        p_xnlength:

    Returns:
    """
    slopeCutoffVertical = 20  # another check

    # Find initial slope:
    if abs(lstThisSegmentCols[0] - lstThisSegmentCols[-1]) < 3:
        m_init = 9999.0
    elif abs(lstThisSegmentRows[0] - lstThisSegmentRows[-1]) < 3:
        m_init = 0.0001
    else:
        m_init = (lstThisSegmentRows[0] - lstThisSegmentRows[-1]) / (
            lstThisSegmentCols[0] - lstThisSegmentCols[-1]
        )

    # Check for zero or infinite slope:
    if m_init == 0:
        m_init = 0.0001
    elif isinf(m_init):
        m_init = 9999.0

    # Find the orthogonal slope:
    m_ortho = -1 / m_init

    xn_steps = [-float(p_xnlength), float(p_xnlength)]  # just the end points

    lst_xy = []
    for r in xn_steps:

        # Make sure it's not too close to vertical:
        # NOTE X-Y vs. Row-Col here:
        if abs(m_ortho) > slopeCutoffVertical:
            tpl_xy = (midPtCol, midPtRow + r)

        else:
            fit_col_ortho = midPtCol + (float(r) / (sqrt(1 + m_ortho**2)))
            tpl_xy = float((midPtCol + (float(r) / (sqrt(1 + m_ortho**2))))), float(
                (m_ortho * (fit_col_ortho - midPtCol) + midPtRow)
            )

        lst_xy.append(tpl_xy)  # A list of two tuple endpts

    return lst_xy


def write_xns_shp(df_coords, epsg, xn_file, p_xngap, xn_type):
    """
    Builds Xns from x-y pairs representing shapely interpolations along a reach

    Args:
        df_coords:
        epsg:
        xn_file:
        bool_isvalley:
        p_xngap:
        logger:

    Returns: list of tuples of lists describing the Xn's along a reach (row, col)
    """
    j = 0

    # slopeCutoffVertical = 20 # just a threshold determining when to call a Xn vertical
    # the final output, a list of tuples of XY coordinate pairs for all Xn's for this reach
    xn_cntr = 0
    lst_xnrowcols = []
    gp_coords = df_coords.groupby("linkno")

    # Create the Xn shapefile for writing:
    test_schema = {
        "geometry": "LineString",
        "properties": {"linkno": "int", "strmord": "int"},
    }

    # logger.info("Building and Writing Cross Section File:")
    with fiona.open(
        xn_file,
        "w",
        driver="ESRI Shapefile",
        crs=f"EPSG:{epsg}",
        schema=test_schema,
    ) as chan_xns:
        for i_linkno, df_linkno in gp_coords:
            i_linkno = int(i_linkno)
            i_order = int(df_linkno.order.iloc[0])
            j += 1

            # NOTE: Define Xn length (p_xnlength) and other parameters relative to stream order
            # Settings for stream channel cross-sections:
            if i_order > 6:
                order_no = 6
            else:
                order_no = i_order
            p_xnlength, p_fitlength = xn_type[str(order_no)].values()

            reach_len = len(df_linkno["x"])

            if reach_len <= p_xngap:
                #                logger.info('Less than!')
                continue  # skip it for now

            # Loop along the reach at the specified intervals:(Xn loop)
            for i in range(p_xngap, reach_len - p_xngap, p_xngap):

                lstThisSegmentRows = []
                lstThisSegmentCols = []

                # if i + paramFitLength > reach_len
                if p_fitlength > i or i + p_fitlength >= reach_len:
                    fitLength = p_xngap
                else:
                    fitLength = p_fitlength

                lstThisSegmentRows.append(df_linkno["y"].iloc[i + fitLength])
                lstThisSegmentRows.append(df_linkno["y"].iloc[i - fitLength])
                lstThisSegmentCols.append(df_linkno["x"].iloc[i + fitLength])
                lstThisSegmentCols.append(df_linkno["x"].iloc[i - fitLength])

                midPtRow = df_linkno["y"].iloc[i]
                midPtCol = df_linkno["x"].iloc[i]

                # Send it the endpts of what you to draw a perpendicular line to:
                lst_xy = build_xns(
                    lstThisSegmentRows,
                    lstThisSegmentCols,
                    midPtCol,
                    midPtRow,
                    p_xnlength,
                )  # returns a list of two endpoints

                xn_cntr = xn_cntr + 1

                # the shapefile geometry use (lon,lat) Requires a list of x-y tuples
                line = {"type": "LineString", "coordinates": lst_xy}
                prop = {"linkno": i_linkno, "strmord": i_order}
                chan_xns.write({"geometry": line, "properties": prop})

    return lst_xnrowcols



def generate(Config, Paths, cell_size):

    coords = get_stream_coords_from_features(Paths.network_poly,
        cell_size,
        Config.preprocess['reach-order']['reach_id'],
        Config.preprocess['reach-order']['order_id'],
        Paths.xn_coordinates
    )

    # channels
    write_xns_shp(coords, Config.spatial_ref['epsg'], Paths.channel_xns, Config.methods['cross_section']['p_xngap'], Config.xn_lengths["channel"])

    # floodplains
    write_xns_shp(coords, Config.spatial_ref['epsg'], Paths.floodplain_xns, Config.methods['cross_section']['p_xngap'], Config.xn_lengths["floodplain"])
