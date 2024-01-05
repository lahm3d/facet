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
from timeit import default_timer as timer

from pathlib import Path

from utils import parse_toml, utils
from utils.batch import generate_processing_batch
from preprocessing import preprocess, cross_sections, network_smoothing
from metrics import channel_cross_section_metrics as channel_metrics
from metrics import channel_curvature_metrics as curvature_metrics
from metrics import flood_inundation_map as fim
from metrics import floodplain_metrics

# debug
# from src.utils import parse_toml, utils
# from src.utils.batch import generate_processing_batch
# from src.preprocessing import preprocess, cross_sections, network_smoothing
# from src.metrics import channel_cross_section_metrics as channel_metrics 
# from src.metrics import channel_curvature_metrics as curvature_metrics 
# from src.metrics import flood_inundation_map as fim
# from src.metrics import floodplain_metrics


# Debug WBT compile issue only on WSL Ubuntu 20.0
# whitebox.download_wbt(linux_musl=True, reset=True)


if __name__ == "__main__":

    config_toml = Path("src/config.toml")
    fpaths_toml = Path("src/utils/filepaths.toml")

    # step 1
    Config = parse_toml.create_config(config_toml)

    # step 2
    hucs = generate_processing_batch(Config.batch_csv)

    for huc in hucs:

        # step 3
        Paths = parse_toml.create_filepaths(fpaths_toml, Config, huc)

        utils.create_folder(Paths)

        # logging
        logger = utils.initialize_logger(Paths.log)

        # Paths_dict = parse_toml.class_to_dict(Paths)
        # for k,v in Paths_dict.items():
        #     print(v)

        # start HUC processing time
        start = timer()

        logger.info(f"Running {huc}...")

        preprocess.run_preprocessing_steps(Config, Paths)

        #### EXPERIMENTAL NETWORK SMOOTHING ####
        # test sample smoothing using 3 refinements
        # untoggle if you want to use smooth version instead of taudem derieved stream network:
        # smooth_network = Paths.network_poly.parent / Paths.network_poly.name.replace('network', 'smooth_network')
        # network_smoothing.apply_chaikins_corner_cutting(Paths.network_poly, smooth_network, refinements=3)
        # Paths.network_poly = smooth_network

        cross_sections.generate(Config, Paths, cell_size=1)

        # 1D Channel Cross-section Metrics
        channel_metrics.derive(
            Paths.channel_xns, Paths.dem, Paths.bank_points,
            Config.methods['cross_section'], Config.spatial_ref['epsg'], logger
            )

        # Channel Curvature Metrics
        curvature_metrics.derive(
            Paths.xn_coordinates, Paths.dem, Paths.bank_pixels, 
            Config.spatial_ref['cell_size'], Config.methods['curvature'], Paths.network_poly,
            Paths.channel_segs, logger
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
            Paths.dem, "ch_wid_tot",
            Paths.channel_segs, Config.xn_lengths["floodplain"],
            logger
        )

        # Experimental HAND method
        floodplain_metrics.hand_method(
            Paths.hand, Paths.channel_segs, Paths.network_poly,
            Paths.network_rast, Paths.flood_extent_layer, Paths.dem,
            Config.xn_lengths["floodplain"], logger
        )

        stop = timer() - start
        print(f"{timer() - start} seconds")
