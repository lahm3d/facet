import sys
from math import isinf, sqrt
from timeit import default_timer as timer

import fiona
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString


def get_stream_coords_from_features(network, xn_gap, min_stream_length, reach_id, order_id, csv_output):
    """
    Generates a large Pandas DataFrame containing a row for every coordinate in
    the stream network. Simplifies the shape of stream reaches and interpolates
    points along each reach at a specified spacing.

    Parameters
    ----------
    network : WindowsPath object of pathlib module
        Path to the vector file defining the stream network, typically a *.shp file,
        previously generated in pre-processing by TauDEM.
    xn_gap : Integer
        Integer defining the spacing between stream network interpolation points,
        in the same linear units as the specified coordinate system.
    min_stream_length : Integer
        Integer defining the minimum length of a stream segment, in the same 
        linear units as the specified coordinate system, in order for the 
        simplify procedure to be carried out.
    reach_id : String
        String defining the name of the attribute associated with reach IDs (LINKNO).
    order_id : String
        String defining the name of the attribute associated with stream order (strmOrder).
    csv_output : WindowsPath object of pathlib module
        Path to a CSV file where the stream coordinates will be written.

    Returns
    -------
    df_final : Pandas DataFrame
        DataFrame containing a row for every coordinate in the stream network with an X-coordinate, Y-coordinate, reach ID, and stream order.

    """
    # Check to see if the coordinates CSV already exists, if so read it and return it instead of re-doing the coordinate generation
    if csv_output.is_file():
        return pd.read_csv(csv_output)
    
    else:
        
        # Enter main coordinate extraction step
        lst_df_final = []
        p_interp_spacing = int( xn_gap )  # 3 # larger numbers would simulate a more smoothed reach
        j = 0  # counter for progress bar (iterating over streamline reaches)
        
        # Open the streamlines shapefile
        with fiona.open(network, "r") as streamlines: 
            tot = len(streamlines)
            
            # Loop over each stream reach simplifying geometry and interpolating coordinates
            for line in streamlines:
                
                j += 1
                line_shply = LineString( line["geometry"]["coordinates"] )
                length = line_shply.length  # units depend on crs

                if length > min_stream_length:  # Skip stream reaches with length less than the specified threshold. NOTE: This value is dependent on CRS!!
                    
                    # Retrieve the reach ID and stream order for the reach
                    i_linkno = line["properties"][reach_id]
                    i_order = line["properties"][order_id]

                    # Smoothing reaches via Shapely, with maximum allowable displacement increasing as a function of stream order:
                    if i_order <= 3:
                        line_shply = line_shply.simplify( 5.0, preserve_topology = False )
                    elif i_order == 4:
                        line_shply = line_shply.simplify( 10.0, preserve_topology = False )
                    elif i_order == 5:
                        line_shply = line_shply.simplify( 20.0, preserve_topology = False )
                    elif i_order >= 6:
                        line_shply = line_shply.simplify( 30.0, preserve_topology = False )

                    length = line_shply.length

                    # Interpolate points from 0 to length of the stream reach at the p_interp_spacing (xn_gap)
                    int_pts = np.arange(0, length, p_interp_spacing)

                    # Generate X/Y coordinates, reach ID, and stream order lists for the interpolated points
                    lst_x = []
                    lst_y = []
                    lst_linkno = []
                    lst_order = []
                    for i in int_pts:
                        i_pt = np.array( line_shply.interpolate( i ) ).item()
                        lst_x.append(i_pt.x)
                        lst_y.append(i_pt.y)
                        lst_linkno.append(i_linkno)
                        lst_order.append(i_order)

                    # Transform the lists into pandas DataFrame
                    df_coords = pd.DataFrame( {"x": lst_x, "y": lst_y,
                                               "LINKNO": lst_linkno, "order": lst_order}
                    )
                    # Remove any potential duplicates caused by interpolation
                    df_coords.drop_duplicates( subset = ["x", "y"], inplace = True )
                    # Append the DataFrame for this reach to the overall DataFrame for all reaches
                    lst_df_final.append( df_coords )
            
            # Concatenate and export the DataFrame of stream network points 
            df_final = pd.concat( lst_df_final )
            df_final.to_csv(csv_output, index = False)

        return df_final


def build_xns(xn_slope_vertical_cutoff, lstThisSegmentRows, lstThisSegmentCols, midPtCol, midPtRow, p_xnlength):
    """
    Generate end points for a single cross-section with a specified length

    Parameters
    ----------
    xn_slope_vertical_cutoff : Intger
        The slope (degrees) cutoff above which a cross section is considered vertical (?).
    lstThisSegmentRows : List
        List of upstream and downstream Y coordinates.
    lstThisSegmentCols : List
        List of upstream and downstream X coordinates.
    midPtCol : float64
        Midpoint X coordinate.
    midPtRow : float64
        Midpoint Y coordinate.
    p_xnlength : Integer
        Length of the cross-section in the same linear units as the specified project coordinate system.

    Returns
    -------
    lst_xy : List
        List of tuples representing either end of the cross section.

    """
    # Calculate the angle (horizontal) from the downstream point to the upstream point (need to confirm order, not sure if reaches go upstream to downstream):
    if abs( lstThisSegmentCols[0] - lstThisSegmentCols[-1] ) < 3:
        m_init = 9999.0
    elif abs(lstThisSegmentRows[0] - lstThisSegmentRows[-1]) < 3:
        m_init = 0.0001
    else:
        m_init = ( lstThisSegmentRows[0] - lstThisSegmentRows[-1] ) / ( 
            lstThisSegmentCols[0] - lstThisSegmentCols[-1] )

    # Check for zero or infinite slope/angle between the upstream and downstream points and replace value if so
    if m_init == 0:
        m_init = 0.0001
    elif isinf( m_init ):
        m_init = 9999.0

    # Calculate the orthogonal slope/angle:
    m_ortho = -1 / m_init

    # Generate a list of the distances to project end points from the  mid point
    xn_steps = [ -float( p_xnlength ), float( p_xnlength ) ]

    # Generate the cross section end points
    lst_xy = []
    for r in xn_steps:

        # Make sure it's not too close to vertical:
        # NOTE X-Y vs. Row-Col here:
        if abs( m_ortho ) > xn_slope_vertical_cutoff:
            # If the angle/slope is too high, just offset the cross-section endpoints horizontally (?)
            tpl_xy = ( midPtCol, midPtRow + r )

        else:
            fit_col_ortho = midPtCol + ( float(r) / ( sqrt( 1 + m_ortho ** 2 ) ) ) 
            tpl_xy = float( ( midPtCol + ( float(r) / ( sqrt(1 + m_ortho ** 2 ) ) ) ) ), float(
                (m_ortho * ( fit_col_ortho - midPtCol ) + midPtRow )
            )

        lst_xy.append(tpl_xy)  # A list of two tuple endpts

    return lst_xy


def write_xns_shp(df_coords, epsg, xn_file, xn_gap, xn_type, xn_slope_vertical_cutoff):
    """
    Constructs cross sections from X-Y coordinate pairs that represent interpolated
    points along a given stream reach.

    Parameters
    ----------
    df_coords : pandas DataFrame
        DESCRIPTION.
    epsg : String
        String specifying the European Petroleum Geospatial Group code defining
        the output horizontal coordinate reference system
    xn_file : WindowPath object of pathlib module
        Path to file where the channel cross sections will be written.
    xn_gap : Integer
        Integer defining the spacing between stream network interpolation points,
        in the same linear units as the specified coordinate system.
    xn_type : Dictionary
        Dictionary mapping specified cross section lengths to stream orders.
    xn_slope_vertical_cutoff : Integer
        DESCRIPTION.

    Returns
    -------
    lst_xnrowcols : TYPE
        list of tuples of lists describing the Xn's along a reach (row, col).

    """
    j = 0 # counter for progress bar (iterating over streamline reaches)

    # xn_slope_vertical_cutoff = 20 # just a threshold determining when to call a Xn vertical
    # the final output, a list of tuples of XY coordinate pairs for all Xn's for this reach
    
    xn_cntr = 0
    lst_xnrowcols = []
    gp_coords = df_coords.groupby("LINKNO") # Group the coordinates DataFrame by reach ID

    # Define a schema for the cross section shapefile:
    test_schema = {
        "geometry": "LineString",
        "properties": {"LINKNO": "int", "strmord": "int"},
    }

    # logger.info("Building and Writing Cross Section File:")
    
    
    # Open up an empty shapefile for writing out the channel cross sections
    with fiona.open( xn_file, "w", driver = "ESRI Shapefile", crs = f"EPSG:{epsg}",
                    schema = test_schema ) as chan_xns:
        
        # Loop over the reach IDs
        for i_linkno, df_linkno in gp_coords:
            
            j += 1
            # Retrieve the reach ID and stream order for the reach
            i_linkno = int(i_linkno)
            i_order = int(df_linkno.order.iloc[0])
     
            # Retrieve cross section length and pf length from the xn_type dictionary (based on reach stream order)
            if i_order > 6:
                order_no = 6
            else:
                order_no = i_order
            p_xnlength, p_fitlength = xn_type[str(order_no)].values()

            # Skip reaches that are shorter than xn_gap and move onto the next reach
            reach_len = len( df_linkno["x"] ) # number of points along reach
            if reach_len <= xn_gap:
                #                logger.info('Less than!')
                continue  # skip it for now

            # Loop across the reach coordinates at the specified interval, (xn_gap)
            for xn_num, i in enumerate( range( xn_gap, reach_len - xn_gap, xn_gap ) ):
                
                # Create empty lists for holding x and y coordinates of stream reach segment line (upstream point, central point, downstream point)
                lstThisSegmentRows = []
                lstThisSegmentCols = []

                # Calculating the fitLength parameter, the distance, in number 
                # of points, to identify upstream and downstream coordinates around
                # the central point
                if p_fitlength > i or i + p_fitlength >= reach_len:
                    # if the pf_len is greater than the distance along the 
                    # reach for the given point, or if the pf_len plus the distance
                    # along the reach for the given point is greater than the
                    # reach length, set fitLenght to the point spacing (xn_gap)
                    fitLength = xn_gap
                else:
                    fitLength = p_fitlength # Else, fitLength =  pf_len
                    
                # Retrieve the X and Y coordinates downstream and upstream of the point of interest, as determined by the fitLength spacing
                lstThisSegmentRows.append( df_linkno[ "y" ].iloc[ i + fitLength ] )
                lstThisSegmentRows.append( df_linkno[ "y" ].iloc[ i - fitLength ] )
                lstThisSegmentCols.append( df_linkno[ "x" ].iloc[ i + fitLength ] )
                lstThisSegmentCols.append( df_linkno[ "x" ].iloc[ i - fitLength ] )
                
                # Retrieve the X and Y coordinates at the point of interest
                midPtRow = df_linkno[ "y" ].iloc[ i ]
                midPtCol = df_linkno[ "x" ].iloc[ i ]
                # Send it the endpts of what you to draw a perpendicular line to:
                
                # Generate cross section
                lst_xy = build_xns(
                    xn_slope_vertical_cutoff,
                    lstThisSegmentRows,
                    lstThisSegmentCols,
                    midPtCol,
                    midPtRow,
                    p_xnlength,
                )  # returns a list of two endpoints

                xn_cntr = xn_cntr + 1

                # Write out the cross sections to a shapefile
                line = {"type": "LineString", "coordinates": lst_xy}
                prop = {"LINKNO": i_linkno, "strmord": i_order}
                chan_xns.write( {"geometry": line, "properties": prop} )

    return lst_xnrowcols


def generate(Config, Paths, logger):
    """
    Generates the network-perpendicular channel and floodplain cross-sections

    Parameters
    ----------
    Config : CreatConfig object of src.utils.parse_tomle module
        Object containing dictionaries specifying processing parameters.
    Paths : CreateFilepaths object of src.utils.parse_tomle module
        Object containing filepaths for the outputs of the FACET workflow.

    Returns
    -------
    None.

    """
    if Config.reuse_xn['reuse_xn_flag'] == False:
        # Generate interpolated coordinates along simplified stream network
        coords = get_stream_coords_from_features(
            Paths.network_poly,
            Config.xn_lengths['xn_gap'],
            Config.xn_lengths['min_length_simplify'],
            Config.preprocess['reach-order']['reach_id'],
            Config.preprocess['reach-order']['order_id'],
            Paths.xn_coordinates
        )
    
        # Generate channel cross sections
        write_xns_shp(
            coords, 
            Config.spatial_ref['epsg'], 
            Paths.channel_xns, 
            Config.xn_lengths['xn_gap'], 
            Config.xn_lengths["channel"],
            Config.xn_lengths['xn_slope_vertical_cutoff'],
        )
    
        # Generate floodplain cross sections
        write_xns_shp(
            coords, 
            Config.spatial_ref['epsg'], 
            Paths.floodplain_xns, 
            Config.xn_lengths['xn_gap'], 
            Config.xn_lengths["floodplain"],
            Config.xn_lengths['xn_slope_vertical_cutoff'],
        )
    else:
        logger.info(f"Skipped generating channel and floodplain cross-sections, reusing cross-sections from {Config.reuse_xn['reuse_xn_data']} version {Config.reuse_xn['reuse_xn_version']}. ")