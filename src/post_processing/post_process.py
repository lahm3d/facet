# -*- coding: utf-8 -*-
"""
Created:            7/1/2019
License:            Creative Commons Attribution 4.0 International (CC BY 4.0)
                    http://creativecommons.org/licenses/by/4.0/
Python version:     Tested on Python 3.7x (x64)


PURPOSE
------------------------------------------------------------------------------
[Floodplain and Channel Evaluation Toolkit: QAQC/Post-processor]

FACET is a standalone Python tool that uses open source modules to map the
floodplain extent and compute stream channel and floodplain geomorphic metrics
such as channel width, streambank height, active floodplain width,
and stream slope from DEMs.

NOTES
------------------------------------------------------------------------------
"""

from functools import reduce
from pathlib import Path
from timeit import default_timer as timer
import multiprocessing as mp
import os
import subprocess
import sys

import geopandas as gpd
import pandas as pd
import numpy as np


def qaqc_R_scripts(r_path, dbf_tuple, logger):
    # get list of R scripts to execute
    dict_rscripts = get_r_scripts_paths()
    bank_pts, chan_segs, fp_xns = dbf_tuple

    for key, value in dict_rscripts.items():
        r_script = value
        r_path, r_script = str(r_path), str(r_script)
        if key == 'install_packages':
            args = [r_path, r_script]
        elif key == 'QAQC_7525_and_9505_bankpoints':
            args = [r_path, r_script, str(bank_pts)]
        elif key == 'QAQC_7525_and_9505_ch_bankpix_2dfp':
            args = [r_path, r_script, str(chan_segs)]
        elif key == 'QAQC_7525_and_9505_fpxns_1d':
            args = [r_path, r_script, str(fp_xns)]
        st = timer()
        # print('printing args   ', args)
        run_r_scripts(args, logger)
        end = round((timer() - st)/60.0, 2)
        logger.info(f'Step 2: Post-processing script {key}. Time elapsed: {end} mins')


def run_r_scripts(args, logger):
    ''' execute R scripts '''
    try:
        subprocess.check_call(args)
    except subprocess.CalledProcessError as e:
        logger.critical(f'Processing {args} script failed. Reason: {e}')
        sys.exit(1)


def get_r_scripts_paths():
    pp_dir = Path(os.path.dirname(os.path.abspath(__file__)))
    install_R = (pp_dir / 'install_packages.r', 'install_packages')
    bank_pts_R = (
        pp_dir / 'QAQC_7525_and_9505_bankpoints.R',
        'QAQC_7525_and_9505_bankpoints'
        )
    twoD_streams_R = (
        pp_dir / 'QAQC_7525_and_9505_ch_bankpix_2dfp.R',
        'QAQC_7525_and_9505_ch_bankpix_2dfp'
        )
    fp_xns_1d_R = (
        pp_dir / 'QAQC_7525_and_9505_fpxns_1d.R',
        'QAQC_7525_and_9505_fpxns_1d'
        )

    raw_list = [install_R, bank_pts_R, twoD_streams_R, fp_xns_1d_R,]
    dict_rscripts = {}

    for i in raw_list:
        f, n = i[0], i[1]
        # check and return file path if true or exit
        if f.is_file():
            dict_rscripts[n] = f
        else:
            print(f"{f}: File does not exist or not found. Please recheck")
            sys.exit(1)

    return dict_rscripts


def clip_lines_inside_poly_mp(arg_tuple):
    '''
    source: https://www.earthdatascience.org/courses/earth-analytics-python/spatial-data-vector-shapefiles/clip-vector-data-in-python-geopandas-shapely/
    function to clip line and polygon data using geopandas
    '''
    # if it encounters a value error because no lines fall inside polygon
    try:
        shp, clip_obj, out_shp = arg_tuple
        # Create a single polygon object for clipping
        poly = clip_obj.geometry.unary_union
        spatial_index = shp.sindex
        # Create a box for the initial intersection
        bbox = poly.bounds
        '''Get a list of id's for each road line that overlaps the bounding box
        and subset the data to just those lines '''
        sidx = list(spatial_index.intersection(bbox))
        shp_sub = shp.iloc[sidx]
        # Clip the data - with these data
        clipped = shp_sub.copy()
        clipped['geometry'] = shp_sub.intersection(poly)
        # Return the clipped layer with no null geometry values
        gdf = clipped[~clipped.geometry.is_empty]
        gdf.to_file(str(out_shp))
    except ValueError:
        return


def flag_FACET_outputs(
        nhd_flowlines, nhd_waterbodies, stream_buffer,
        huc10_mask, twoD_streams, bank_pts,
        fp_xns, out_dir, logger):

    '''
    QAQC flagging of bank points and stream xns that are either
    outside of NHD flowline buffer OR inside NHD waterbodies

    inputs:
    1. NHD flowlines
    2. NHD Waterbodies
    3. HUC10 boundary or mask
    4. FACET: 2D stream network
    5. FACET: bank points
    6. FACET: floodplain x-sections

    outputs:
    New 2d-streams, fp xns and bank pts shps where features that
    fall ouside the NHD buffer or water bodies are flagged
    '''
    if not out_dir.is_dir():
        os.mkdir(out_dir)

    st = timer()
    logger.info('Step 1: Post-processing...')

    # NHD derived hydrographic files
    gdf_flowlines = gpd.read_file(str(nhd_flowlines))
    gdf_wb = gpd.read_file(str(nhd_waterbodies))
    gdf_mask = gpd.read_file(str(huc10_mask))
    # FACET derived files
    gdf_2D_streams = gpd.read_file(str(twoD_streams))
    gdf_bank_pts = gpd.read_file(str(bank_pts))
    gdf_fp_xns = gpd.read_file(str(fp_xns))

    # tmp files
    tmp_outside_buffer = out_dir / 'tmp_outside_buffer.shp'
    tmp_water_bodies = out_dir / 'tmp_water_bodies.shp'
    tmp_streams_outside_nhd_buff = out_dir / 'tmp_streams_outside_nhd_buff.shp'
    tmp_streams_inside_waterbodies = out_dir / 'tmp_streams_inside_waterbodies.shp'
    out_bank_pts = out_dir / 'bankpoints_1D_metrics.shp'
    out_fp_xns = out_dir / 'floodplain_xns_1D_metrics.shp'
    out_2D_streams = out_dir / 'channel_floodplain_2D_metrics.shp'

    # 0 - generate unique IDs to perform joins
    uid = 'unique_id'
    gdf_2D_streams[uid] = np.arange(gdf_2D_streams.shape[0])
    gdf_bank_pts[uid] = np.arange(gdf_bank_pts.shape[0])
    gdf_fp_xns[uid] = np.arange(gdf_fp_xns.shape[0])

    # 1 - Universe of buffered NHD streams
    '''create inverse polygon of NHD buffered streams &
        clip streams by outside buffer'''
    logger.info('starting: create outside buffer mask & streams outside buffer')
    nhd_buffer = gdf_flowlines.copy()
    nhd_buffer['geometry'] = nhd_buffer.buffer(stream_buffer)
    nhd_buffer['id'] = 1
    nhd_buffer = nhd_buffer.dissolve(by='id')

    outside_buffer = gpd.overlay(gdf_mask, nhd_buffer, how='symmetric_difference')
    outside_buffer.to_file(str(tmp_outside_buffer))

    # 2 - query NHD water bodies & clip streams inside waterbodies
    ''' FType codes for lake/pond is 390 and reservoir 436 '''
    water_bodies = gdf_wb.loc[gdf_wb['FType'].isin([390, 436])]
    water_bodies = gpd.overlay(water_bodies, gdf_mask, how='intersection')
    water_bodies.to_file(str(tmp_water_bodies))

    # 3 - clip streams inside outside_buffer and water bodies
    clip_list = [
        (gdf_2D_streams, outside_buffer, tmp_streams_outside_nhd_buff),
        (gdf_2D_streams, water_bodies, tmp_streams_inside_waterbodies)
        ]
    pool = mp.Pool(processes=2)
    pool.map(clip_lines_inside_poly_mp, [(item) for item in clip_list])
    pool.close()

    # check to see if files exist
    polygon_list = []
    if tmp_streams_outside_nhd_buff.is_file():
        streams_outside_buff = gpd.read_file(str(tmp_streams_outside_nhd_buff))
        polygon_list.append((streams_outside_buff.copy(), 'NHDFlag'))

    if tmp_streams_inside_waterbodies.is_file():
        streams_inside_wb = gpd.read_file(str(tmp_streams_inside_waterbodies))
        polygon_list.append((streams_inside_wb.copy(), 'WBDFlag'))


    # 5 - Flag bank points inside water bodies & within outside-NHD-flowline buffer
    logger.info('Flagging bank points:')
    bank_pts_gdfs = []
    for poly in [(outside_buffer.copy(), 'NHDFlag'), (water_bodies.copy(), 'WBDFlag')]:
        poly_gdf, col_name = poly # split tuples
        logger.info(f'\tstart intersection with ch xns: {col_name}')
        bank_pts = gdf_bank_pts.copy()
        gdf = gpd.sjoin(bank_pts, poly_gdf, how='inner', op='within')
        list_of_xn_nums = list(set(gdf.xn_num))
        gdf = bank_pts.loc[bank_pts['xn_num'].isin(list_of_xn_nums)]
        gdf = gdf.drop_duplicates(subset=[uid], keep='first')
        gdf[col_name] = 1 # flag
        gdf = gdf[[uid, col_name]] # clean columns
        logger.info(f'\tappend: {col_name}')
        bank_pts_gdfs.append(gdf)

    # merge data frames
    # extract gdf from list with specific columns only
    gdf1, gdf2 = bank_pts_gdfs[0][[uid, 'NHDFlag']], bank_pts_gdfs[1][[uid, 'WBDFlag']]
    # sequential merges
    merge_1 = gdf_bank_pts.merge(gdf1, how='left', on=uid)
    final_bank_pts = merge_1.merge(gdf2, how='left', on=uid)
    logger.info('Bank points successfully flagged\n')

    # 6 - Flag fp xsections that intersect streams inside water bodies & outside-NHD-flowline buffer
    logger.info('Flagging Floodplain Cross-sections')
    fp_xns_gdfs = []
    for poly in polygon_list:
        poly_gdf, col_name = poly # split tuples
        logger.info(f'\tstart overlay: {col_name}')
        poly_gdf['strm_uid'] = np.arange(poly_gdf.shape[0])
        poly_gdf = poly_gdf[['strm_uid', 'geometry']]
        gdf = gpd.sjoin(gdf_fp_xns.copy(), poly_gdf, how='inner', op='intersects')
        gdf = gdf.drop_duplicates(subset=[uid], keep='first')
        gdf[col_name] = 1 # flag
        logger.info(f'\tappend: {col_name}')
        fp_xns_gdfs.append(gdf)

    # merge data frames
    # extract gdf from list with specific columns only
    if len(fp_xns_gdfs) == 1:
        gdf = fp_xns_gdfs[0][[uid, polygon_list[0][1]]] # get flag from polygon_list
        final_fp_xns = gdf_fp_xns.merge(gdf, how='left', on=uid)
    else:
        gdf1, gdf2 = fp_xns_gdfs[0][[uid, 'NHDFlag']], fp_xns_gdfs[1][[uid, 'WBDFlag']]
        # sequential merges
        merge_1 = gdf_fp_xns.merge(gdf1, how='left', on=uid)
        final_fp_xns = merge_1.merge(gdf2, how='left', on=uid)
    logger.info('Floodplain Cross-sections successfully flagged\n')

    # 7 - Flag streams inside water bodies and streams outside the buffer
    logger.info('Flagging 2D streams')
    flagged_2D_strm_gdfs = []
    for poly in polygon_list:
        poly_gdf, col_name = poly # split tuples
        logger.info(f'\tstart flagging: {col_name}')
        poly_gdf[col_name] = 1 # flag
        logger.info(f'\tappend: {col_name}')
        poly_gdf = poly_gdf[[uid, col_name]] # clean columns
        flagged_2D_strm_gdfs.append(poly_gdf)

    # merge data frames
    flagged_2D_streams = reduce(
        lambda left, right: pd.merge(left, right, how='outer', on=uid),
        flagged_2D_strm_gdfs
        )
    final_2D_streams = gdf_2D_streams.merge(flagged_2D_streams, how='left', on=uid)
    logger.info('2D streams successfully flagged\n')


    # 8 - write out flagged fp-xns and bank pts
    final_bank_pts.to_file(str(out_bank_pts))
    final_fp_xns.to_file(str(out_fp_xns))
    final_2D_streams.to_file(str(out_2D_streams))

    end = round((timer() - st)/60.0, 2)
    logger.info(f'Step 1: Post-processing complete. Time elapsed: {end} mins')

    out_tuple = (
        out_dir / 'bankpoints_1D_metrics.dbf',
        out_dir / 'channel_floodplain_2D_metrics.dbf',
        out_dir /  'floodplain_xns_1D_metrics.dbf'
        )

    return out_tuple


def rename_shp_fields(huc_dir, logger):
    """ rename fields based on file dictionary to create final outputs """

    # list of SHPs to rename
    post_process_shps = (
        huc_dir / 'bankpoints_1D_metrics.shp',
        huc_dir / 'channel_floodplain_2D_metrics.shp',
        huc_dir / 'floodplain_xns_1D_metrics.shp',
        )

    # dict of files with old and new field name
    fields_dict = {
        'bankpoints_1D_metrics': {
            'bank_hght': 'bankht_1d',
            'bf_area': 'chan_area',
        },
        'floodplain_xns_1D_metrics': {
            'strmord': 'strmorder',
            'xn_id_1dfp': 'fpxn_1d',
            'totwid_1df': 'fpwid_1d',
            'mindep_1df': 'mind_1d',
            'maxdep_1df': 'maxd_1d',
            'rngdep_1df': 'rngd_1d',
            'meandep_1d': 'meand_1d',
            'stddep_1df': 'stdd_1d',
            'sumdep_1df': 'sumd_1d',
            'minele_1df': 'mine_1d',
            'maxele_1df': 'maxe_1d',
            'rngele_1df': 'rnge_1d',
            'meanele_1d': 'meane_1d',
            'stdele_1df': 'stde_1d',
            'sumele_1df': 'sume_1d',
            },
        'channel_floodplain_2D_metrics': {
            'ch_wid_tot': 'chnwid_px',
            'ch_wid_1': 'chnwid1_px',
            'ch_wid_2': 'chnwid2_px',
            'order': 'strmorder',
            'fp_width_2': 'fpwid_2dc',
            'fp_range_2': 'rngd_2dc',
            'fp_min_2d': 'mind_2dc',
            'fp_max_2d': 'maxd_2dc',
            'fp_std_2d': 'stdd_2dc',
            'fp_rug_2d': 'rug_2dc',
            'bnk_ht3': 'bankht_2dh',
            'chn_shp3': 'chnshp_2dh',
            'chn_wid3': 'chnwid_2dh',
            'fp3_min': 'mind_2dh',
            'fp3_max': 'maxd_2dh',
            'fp3_std': 'stdd_2dh',
            'fp3_wid': 'fpwid_2dh',
            'fp3_rug': 'rug_2dh',
            'fp3_rng': 'rngd_2dh',
            'fp3_min_e': 'mine_2dh',
            'fp3_max_e': 'maxe_2dh',
            'fp3_std_e': 'stde_2dh',
            'fp3_rng_e': 'rnge_2dh',
            },
    }

    # rename columns in files based on the dict
    for i in post_process_shps:
        gdf = gpd.read_file(str(i))
        bname = i.stem
        try:
            if set(fields_dict[bname].keys()).issubset(gdf.columns):
                old_cols = set(gdf.columns)
                gdf.rename(columns=fields_dict[bname], inplace=True, errors="raise")
                diff_cols = set(gdf.columns) - old_cols # diff should be new field names
                new_cols = set(fields_dict[bname].values()) # get new cols as set
                # additional check to see if new fields are correct
                if diff_cols == new_cols:
                    gdf.to_file(str(i))
                    logger.info(f'Fields successfully renamed {bname}: {diff_cols}')
        except KeyError as e:
            logger.warning(f'{i} encountered error: {e}')
