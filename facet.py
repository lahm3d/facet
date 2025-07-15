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

------------------------------------------------------------------------------
"""
import argparse
import time
from pathlib import Path

from src.utils import parse_toml, utils
from src.utils.batch import generate_processing_batch
from src.preprocessing import preprocess, cross_sections, network_smoothing
from src.metrics import channel_cross_section_metrics as channel_metrics
from src.metrics import channel_curvature_metrics as curvature_metrics
from src.metrics import flood_inundation_map as fim
from src.metrics import floodplain_metrics

from src.postprocessing import spatial_qc as qc


# Debug WBT compile issue only on WSL Ubuntu 20.0
# whitebox.download_wbt(linux_musl=True, reset=True)

if __name__ == "__main__":
    
    # Command line argument parsing for easier command line implementation
    parser = argparse.ArgumentParser()
    parser.add_argument("--config_toml", help = "filepath to configuration toml file, relative to root FACET directory, i.e., src/config_test.toml")
    parser.add_argument("--fpaths_toml", help = "filepath to filepaths toml file, relative to root FACET directory, i.e. src/filepaths.toml")
    args = parser.parse_args()
    
    config_toml = Path(args.config_toml)
    fpaths_toml = Path(args.fpaths_toml)
    
    # config_toml = Path("src/config_test.toml")
    # fpaths_toml = Path("src/utils/filepaths.toml")

    # step 1
    Config = parse_toml.create_config(config_toml)

    # step 2
    if ( Config.hucs['batch_csv'] != "None" ) & ( Config.hucs['huc'] != "None" ):
        raise ValueError(f'Both a CSV of HUCs and an individual HUC are specified in the .toml, choose one or the other')
    elif ( Config.hucs['batch_csv'] == "None" ) & ( Config.hucs['huc'] == "None" ):
        raise ValueError(f'Both a CSV of HUCs and an individual HUC are set to "None" in the .toml, one these needs to have a valid value')
    elif ( Config.hucs['batch_csv'] != "None" ):
        hucs = generate_processing_batch(Config.batch_csv)
    elif ( Config.hucs['huc'] != "None" ):
        hucs = [ Config.hucs['huc'] ] # create a list containing the single HUC code so that it is iterable and the below for-loop does not break

    for huc in hucs:
        # step 3
        Paths = parse_toml.create_filepaths(fpaths_toml, Config, huc)

        utils.create_folder(Paths)

        # logging
        logger = utils.initialize_logger(Paths.log)

        # log input parameters
        Paths_dict = parse_toml.class_to_dict(Paths)
        for k,v in Paths_dict.items():
            logger.debug(f"{k}: {v}")

        # start HUC processing time
        start = time.time()

        logger.info(f"Running {huc}...")

        preprocess.run_preprocessing_steps(Config, Paths, logger)

        #### EXPERIMENTAL NETWORK SMOOTHING ####
        # test sample smoothing using 3 refinements
        # untoggle if you want to use smooth version instead of taudem derieved stream network:
        # smooth_network = Paths.network_poly.parent / Paths.network_poly.name.replace('network', 'smooth_network')
        # network_smoothing.apply_chaikins_corner_cutting(Paths.network_poly, smooth_network, refinements=3)
        # Paths.network_poly = smooth_network

        # Generate cross-sections
        cross_sections.generate(Config, Paths, logger)

        # 1D Channel Cross-section Metrics
        channel_metrics.derive(
            Config.methods['cross_section']['cell_size'],
            Paths.elevation_profiles,
            Paths.channel_xns,
            Paths.dem,
            Paths.bank_points,
            Config.methods['cross_section'],
            Config.spatial_ref['epsg'],
            logger
            )

        # Channel Curvature Metrics
        curvature_metrics.derive(
            Paths.xn_coordinates,
            Paths.dem,
            Paths.bank_pixels,
            Config.spatial_ref['cell_size'],
            Config.methods['curvature'],
            Paths.network_poly,
            Paths.channel_segs,
            logger
            )

        # delineate flood inundation layer
        reach_id = Config.preprocess['reach-order']['reach_id']
        min_da, max_da = Config.methods['flood_thresholds'].values()
        fim.delineate(
            Paths.hand,
            Paths.sub_watersheds_poly,
            Config.preprocess['reach-order']['reach_id'],
            Paths.flood_extent_layer,
            Paths.flood_height_thresholds,
            min_da,
            max_da,
            logger
            )

        # 1D Floodplain Cross-section Metrics
        floodplain_metrics.derive(
            Paths.floodplain_xns, Paths.flood_extent_layer,
            Paths.dem, "CURVE_WD",
            Paths.channel_segs, Config.xn_lengths["floodplain"],
            logger
        )

        # # Experimental HAND method
        # floodplain_metrics.hand_method(
        #     Paths.hand, Paths.channel_segs, Paths.network_poly,
        #     Paths.network_rast, Paths.flood_extent_layer, Paths.dem,
        #     Config.xn_lengths["floodplain"], logger
        # )

        # Quality Checks against NHD
        flowline_mask = qc.create_flowline_qc_mask(
            Paths.flowlines, 
            Config.postprocess['stream-buffer'],
            Paths.watershed
        )

        waterbody_mask = qc.create_waterbody_qc_mask(
            Config.ancillary['nhd_wbds'],
            [390, 436], 
            Paths.watershed
            )

        # flag bankpoints:
        bank_points_qc = utils.vector_to_geodataframe(Paths.bank_points)
        bank_points_qc = qc.flag_features_by_qc_mask(
            bank_points_qc, flowline_mask, "NHD_Flag", "xn_num"
            )
        bank_points_qc = qc.flag_features_by_qc_mask(
            bank_points_qc, waterbody_mask, "WBD_Flag", "xn_num", output=Paths.bank_points
            )

        # flag channel segs:
        channel_segs_qc = utils.vector_to_geodataframe(Paths.channel_segs)

        channel_segs_qc = qc.flag_features_by_qc_mask(
            channel_segs_qc, flowline_mask, "NHD_Flag"
            )
        channel_segs_qc = qc.flag_features_by_qc_mask(
            channel_segs_qc, waterbody_mask, "WBD_Flag", output=Paths.channel_segs
            )

        # flag floodplain xns:
        floodplain_xns_qc = utils.vector_to_geodataframe(Paths.floodplain_xns)

        floodplain_xns_qc = qc.flag_features_by_qc_mask(
            floodplain_xns_qc, flowline_mask, "NHD_Flag"
            )
        floodplain_xns_qc = qc.flag_features_by_qc_mask(
            floodplain_xns_qc, waterbody_mask, "WBD_Flag", output=Paths.floodplain_xns
            )

        logger.info(f"Total run time: {round((time.time() - start) / 60, 2)} mins")
