# -*- coding: utf-8 -*-
"""
Created:            6/7/2019
License:            Creative Commons Attribution 4.0 International (CC BY 4.0)
                    http://creativecommons.org/licenses/by/4.0/
Python version:     Tested on Python 3.7x (x64)


PURPOSE
------------------------------------------------------------------------------
[Floodplain and Channel Evaluation Toolkit]

small module to change FACET parameters/config

NOTES
------------------------------------------------------------------------------
"""

from pathlib import Path
import os
import sys
import configparser
import pandas as pd


def get_config_path():
    ''' this retrieves config.ini provided with the repo '''
    config_dir = Path(os.path.dirname(os.path.abspath(__file__))) # get script directory
    file = config_dir / "config.ini" # construct config file path

    # check and return file path if true or exit
    if file.exists() is True:
        return file
    else:
        print("File does not exist or not found. Please recheck!")
        sys.exit(1)


def validate_breach_methods(cf_dict):
    list_of_methods = [
        cf_dict['pre-condition dem & fast-breach'],
        cf_dict['fast-breach']
        ]
    method_sum = sum(map(int, list_of_methods))

    ''' default behavior if multiple methods are turned on
    then 'pre-condition dem & fast-breach' is selected' '''

    if method_sum == 1:
        return
    else:
        print('Multiple methods selected default set to pre-condition dem & fast-breach')
        cf_dict['pre-condition dem & fast-breach'] = True


def validate_resolution(cf_dict):
    """ sets resample resolution as xnptdist """
    # get pixel size
    pixelSize = int(cf_dict['resample resolution'])

    # ensure resample resolution is same as xnptdist
    cf_dict['xnptdist'] = pixelSize

    # update resample resolution formating
    cf_dict['resample resolution'] = (pixelSize, pixelSize)


def validate_skip_huc10_list(cf_dict):
    if not cf_dict['skip_list']:
        cf_dict['skip_list'] = ('No values provided', 'No validation test')
    else:
        string_2_list = cf_dict['skip_list'][0].split(',')
        cf_dict['skip_list'] = (string_2_list, 'pass')


def validate_config_params(config_file):
    ''' Config file is passed as input and each item is validated '''

    # check if file exists and read it
    if Path(config_file).is_file():
        Config = configparser.ConfigParser()
        Config.read(config_file)
    else:
        print(f'{config_file} file not found!')

    # capture output params
    cf_dict = {}

        # data type lists
    int_params = [
        'taudem cores', 'stream buffer', 'resample resolution',
        'p_xngap', 'i_step', 'max_buff',]

    float_params = [
        'parm_ivert', 'parm_ratiothresh', 'parm_slpthresh']

    bool_params = [
        'taudem',
        'wt_grid',
        'pre-condition dem & fast-breach',
        'fast-breach',
        'use_wavelet_curvature_method',
        'post process',
        'clean',
        'archive']

    path_params = [
        'data_dir', 'ancillary dir', 'physio cbw',
        'physio drb', 'census roads', 'census rails',
        'nhd_dir', 'r exe path']

    cf_dict = {}
    # validate correct passing correct params
    for section in Config.sections():
        for param in Config[section]:
            value = Config.get(section, param)
            # fix default validation parameters
            if param in bool_params:
                value = Config.getboolean(section, param)
            elif param in int_params:
                value = Config.getint(section, param)
            elif param in float_params:
                value = Config.getfloat(section, param)
            elif param in path_params:
                value = Path(Config.get(section, param))
            # set key and value
            cf_dict[param] = value

    # validate only one breach method is selected
    validate_breach_methods(cf_dict)

    # update resample resolution and xnptdist
    validate_resolution(cf_dict)

    # check if file passed for skip huclist or string
    try:
        df = pd.read_csv(cf_dict['skip_list'], header=None)
        df[0] = '0' + df[0].astype(str)
        cf_dict['skip_list'] = df[0].to_list()
    except FileNotFoundError:
        cf_dict['skip_list'] = cf_dict['skip_list'].split(',')

    return cf_dict
