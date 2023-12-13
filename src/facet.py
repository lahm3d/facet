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
from preprocessing import preprocess, cross_sections

from src.utils import parse_toml, utils
from src.utils.batch import generate_processing_batch
from src.preprocessing import preprocess, cross_sections

# import pandas as pd
# from configure import config
# from post_processing import post_process
# import preprocessing
# import utils.utils as utils

# from metrics import (
#     channel_cross_section_metrics,
#     floodplain_metrics,
#     channel_curvature_metrics,
# )

if __name__ == "__main__":


    config_toml = Path("src/config_test.toml")
    fpaths_toml = Path("src/utils/filepaths.toml")

    # step 1
    Config = parse_toml.create_config(config_toml)

    # step 2
    hucs = generate_processing_batch(Config.batch_csv)

    for huc in hucs:

        # step 3
        Paths = parse_toml.create_filepaths(fpaths_toml, Config, huc)

        utils.create_folder(Paths)

        # Paths_dict = parse_toml.class_to_dict(Paths)
        # for k,v in Paths_dict.items():
        #     print(v)
        
        # start HUC processing time
        start = timer()

        # set-up logging
        # log_file, time_stamp = utils.setup_logging(huc, Paths.parent)
        # logger = utils.initialize_logger(log_file)
        # logger.info(f"Running {HUCID}. Start time stamp: {time_stamp}")

        preprocess.run_preprocessing_steps(Config, Paths)

        cross_sections.generate(Config, Paths, cell_size=1)

        stop = timer() - start
        print(f"{timer() - start} seconds")

