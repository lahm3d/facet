# -*- coding: utf-8 -*-
"""
Created:            6/7/2019
License:            Creative Commons Attribution 4.0 International (CC BY 4.0)
                    http://creativecommons.org/licenses/by/4.0/
Python version:     Tested on Python 3.7x (x64)


PURPOSE
------------------------------------------------------------------------------
[Floodplain and Channel Evaluation Toolkit]

FACET is a standalone Python tool that uses open source modules to map the
floodplain extent and compute stream channel and floodplain geomorphic metrics
such as channel width, streambank height, active floodplain width,
and stream slope from DEMs.

NOTES
------------------------------------------------------------------------------
"""

from pathlib import Path
from timeit import default_timer as timer
import sys

import pandas as pd
from setup import config
import funcs
from post_processing import post_process

if __name__ == "__main__":
    print("\n<<< Start >>>\r\n")

    # read in config file & generate PARAMS dict
    CONFIG_FILE = config.get_config_path()
    PARAMS = config.validate_config_params(CONFIG_FILE)

    # generate list of huc folders in data dir
    HUC_LST = []
    for p in PARAMS["data_dir"].iterdir():
        if p.is_dir():
            HUC_LST.append(p)

    # loop through the HUC folders
    for huc_dir in HUC_LST:
        HUCID = str(huc_dir).split("\\")[-1]

        # start HUC processing time
        start_time_i = timer()

        # ========== Initialize Logging ========== #
        # set-up logging
        log_file, time_stamp = funcs.setup_logging(HUCID, huc_dir)

        # initialize logging
        logger = funcs.initialize_logger(log_file)

        # first log entry
        logger.info(f"Running {HUCID}. Start time stamp: {time_stamp}")

        # ========== Set-up Ancillary Datasets ========== #
        HUC04 = HUCID[0:4]  # HUC 04 ID

        # Reproject Ancillary Data
        funcs.reproject_ancillary_data(HUC04, PARAMS, logger)

        # skip HUC if in skip list
        if HUCID in PARAMS["skip_list"]:
            logger.warning(f"{HUCID} was excluded based on the skip list")
            continue

        # ========== Start Processing HUC ========== #
        """
        1. Construct file paths and folders each HUC
        2. Check *_DEM.tif and *_mask.shp in each HUC and reproject
        3. Pre-processing steps
             i. Pre-condition DEM for breach
            ii. TauDEM Pre-processing steps
        4. Core FACET Algorithms
        5. Post-processing
        """

        # File paths and Folders for each HUC
        str_nhdhr_huc10 = huc_dir / f"{HUCID}_dem_nhdhires.shp"  # NHD hi-res flowlines
        str_dem_proj = huc_dir / f"{HUCID}_dem_proj.tif"  # projected DEM

        # pre-processing file paths
        str_breached_dem_path = huc_dir / f"{HUCID}_breach.tif"  # breached DEM
        str_hand_path = huc_dir / f"{HUCID}_hand.tif"
        str_net_path = huc_dir / f"{HUCID}_network.shp"
        str_raster_net_path = huc_dir / f"{HUCID}_network.tif"
        str_sheds_path = huc_dir / f"{HUCID}_breach_w_diss_physio.shp"

        # core FACET file paths
        str_csv_path = huc_dir / f"{HUCID}.csv"
        str_chxns_path = huc_dir / f"{HUCID}_channel_xns.shp"
        str_bankpts_path = huc_dir / "bankpoints_1D_metrics.shp"
        str_chanmet_segs = huc_dir / "channel_floodplain_2D_metrics.shp"
        str_bankpixels_path = huc_dir / f"{HUCID}_bankpixels.tif"
        str_fpxns_path = huc_dir / "floodplain_xns_1D_metrics.shp"
        str_fim_path = huc_dir / f"{HUCID}_floodplain.tif"
        str_fim_csv = huc_dir / f"{HUCID}_floodplain_hand_height.csv"

        # FACET post-processing QAQC
        post_process_dir = huc_dir / "post_processing"

        # ========== FACET Pre-processing Steps ========== #
        # raw dem file path
        str_dem = huc_dir / f"{HUCID}_dem.tif"

        # Check if DEM is exists
        if str_dem.is_file():
            logger.info(f"{HUCID}: DEM exists.")
        else:
            logger.critical(f"{HUCID}: DEM NOT found. This HUC will be skipped.")
            continue

        # Check if mask is exists
        try:
            huc10_mask = list(huc_dir.rglob("*_mask.shp"))[0]  # raw HUC10_mask
            # reproject the mask
            huc10_mask_proj = funcs.reproject_vector_layer(
                huc10_mask, PARAMS["crs"], logger
            )
            logger.info(f"{HUCID}: Mask SHP exists.")
        except IndexError:
            logger.critical(f"{HUCID}: Mask SHP NOT found. This HUC will be skipped.")
            continue

        # Project dem raster
        if not str_dem_proj.is_file():
            str_dem_path_proj = funcs.reproject_grid_layer(
                str_dem,
                PARAMS["crs"],
                str_dem_proj,
                PARAMS["resample resolution"],
                logger,
            )
        else:
            logger.info(
                f"{str_dem_proj} already exists. DEM re-projection and resampling skipped."
            )
            str_dem_path_proj = str_dem_proj

        # Clip NHD plus HR streams by HUC:
        funcs.clip_nhdhr_using_grid(
            PARAMS["streams prj"],
            str_nhdhr_huc10,
            str_dem_path_proj,
            PARAMS["crs"],
            logger,
            huc10_mask_proj,
        )

        # pre-breach conditioning + breaching
        if PARAMS["pre-condition dem & fast-breach"]:
            # default breach
            dem_merge = funcs.pre_breach_DEM_conditioning(
                huc_dir,
                HUCID,
                str_dem_proj,
                str_nhdhr_huc10,
                PARAMS["census roads prj"],
                PARAMS["census rails prj"],
                huc10_mask_proj,
                logger,
            )
            funcs.breach_dem(dem_merge, str_breached_dem_path, logger)
        elif PARAMS["fast-breach"]:
            funcs.breach_dem(str_dem_proj, str_breached_dem_path, logger)

        # TauDEM preprocessing steps
        funcs.preprocess_dem(
            huc_dir,
            str_nhdhr_huc10,
            PARAMS["crs"],
            PARAMS["wt_grid"],
            PARAMS["taudem"],
            PARAMS["physio prj"],
            HUCID,
            str_breached_dem_path,
            PARAMS["taudem cores"],
            logger,
        )

        # ========== Core FACET Algorithms ========== #
        st_core_facet_time = timer()

        # Convert vector streamlines to raster with pixel streamline values matching linkno:
        funcs.rasterize_gdf(str_net_path, str_hand_path, str_raster_net_path)

        # << GET CELL SIZE >>
        cell_size = int(funcs.get_cell_size(str_dem_proj))  # range functions need int?

        # << BUILD STREAMLINES COORDINATES >>
        logger.info("Generating the stream network coordinates from the csv file...")
        df_coords, streamlines_crs = funcs.get_stream_coords_from_features(
            str_net_path, cell_size, PARAMS["reach_id"], PARAMS["order_id"], logger
        )
        df_coords.to_csv(str_csv_path)
        logger.info("Reading the stream network coordinates from the csv file...")
        df_coords = pd.read_csv(
            str_csv_path,
        )

        # ========== << CROSS SECTION ANALYSES >> ==========
        # << CREATE Xn SHAPEFILES >>
        # Channel:
        funcs.write_xns_shp(
            df_coords, streamlines_crs, str_chxns_path, False, PARAMS["p_xngap"], logger
        )
        # Floodplain:
        funcs.write_xns_shp(
            df_coords, streamlines_crs, str_fpxns_path, True, int(30), logger
        )

        # << INTERPOLATE ELEVATION ALONG Xns >>
        df_xn_elev = funcs.read_xns_shp_and_get_dem_window(
            str_chxns_path, str_dem_proj, logger
        )

        # Calculate channel metrics and write bank point shapefile
        funcs.chanmetrics_bankpts(
            df_xn_elev,
            str_chxns_path,
            str_dem_proj,
            str_bankpts_path,
            PARAMS["parm_ivert"],
            PARAMS["xnptdist"],
            PARAMS["parm_ratiothresh"],
            PARAMS["parm_slpthresh"],
            logger,
            PARAMS["crs"],
        )

        # ========== << BANK PIXELS AND WIDTH FROM CURVATURE >> ==========
        funcs.bankpixels_from_curvature_window(
            df_coords,
            str_dem_proj,
            str_bankpixels_path,
            cell_size,
            PARAMS["use_wavelet_curvature_method"],
            logger,
        )

        funcs.channel_width_from_bank_pixels(
            df_coords,
            str_net_path,
            str_bankpixels_path,
            PARAMS["reach_id"],
            PARAMS["i_step"],
            PARAMS["max_buff"],
            str_chanmet_segs,
            logger,
        )

        # ========== << DELINEATE FIM >> ==========
        funcs.fim_hand_poly(
            str_hand_path,
            str_sheds_path,
            PARAMS["reach_id"],
            str_fim_path,
            str_fim_csv,
            logger,
        )

        # ========== << FLOODPLAIN METRICS >> ==========
        # 1D approach:
        funcs.read_fp_xns_shp_and_get_1D_fp_metrics(
            str_fpxns_path, str_fim_path, str_dem_proj, logger
        )

        # 2D approach:
        funcs.fp_metrics_chsegs(str_fim_path, "ch_wid_tot", str_chanmet_segs, logger)

        # ========== << HAND CHARACTERISTICS >> ==========

        # Calculate channel and FP metrics through analysis of the HAND grid, separating the
        # in-channel pixels from the FP pixels.

        funcs.hand_analysis_chsegs(
            str_hand_path,
            str_chanmet_segs,
            str_raster_net_path,
            str_fim_path,
            str_dem_proj,
            logger,
        )

        # capture core run time for FACET algorithms
        end_core_facet_time = round((timer() - st_core_facet_time) / 60.0, 2)
        logger.info(
            f"Core FACET algorithms completed. Time elapsed: {end_core_facet_time} mins"
        )

        # ========== Post-Processing  ==========
        # rename columns in shapefile
        post_process.rename_shp_fields(huc_dir, logger)

        if PARAMS["post process"]:
            dbf_tuple = post_process.flag_FACET_outputs(
                str_nhdhr_huc10,
                PARAMS["waterbody prj"],
                PARAMS["stream buffer"],
                huc10_mask_proj,
                str_chanmet_segs,
                str_bankpts_path,
                str_fpxns_path,
                post_process_dir,
                logger,
            )

            post_process.qaqc_R_scripts(PARAMS["r exe path"], dbf_tuple, logger)

        # ========== Clean-up ==========
        if PARAMS["clean"]:
            funcs.clean_tmp_files(huc_dir)

        if PARAMS["archive"]:
            funcs.archive_huc(huc_dir)

        end_time = round((timer() - start_time_i) / 60.0, 2)
        logger.info(f"\nAll steps complete for {HUCID}:  {end_time} mins\r\n")
