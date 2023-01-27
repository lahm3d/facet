# -*- coding: utf-8 -*-
"""
Created:            7/1/2019
License:            Creative Commons Attribution 4.0 International (CC BY 4.0)
                    http://creativecommons.org/licenses/by/4.0/
Python version:     Tested on Python 3.7x (x64)


PURPOSE
------------------------------------------------------------------------------
[Floodplain and Channel Evaluation Toolkit: Setup Script]

FACET is a standalone Python tool that uses open source modules to map the
floodplain extent and compute stream channel and floodplain geomorphic metrics
such as channel width, streambank height, active floodplain width,
and stream slope from DEMs.

------------------------------------------------------------------------------
"""

from pathlib import Path
from timeit import default_timer as timer
import argparse
import configparser
import json
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import zipfile
import requests

import config

def json_to_dict(project_dir):
    """
    Reads science base json listing for FACET Ancillary Data
    and returns a dictionary of urls of ancillary datasets

    The datasets parsed and organized based on their type into
    a new dictionary
    """
    json_urls = [
        "https://www.sciencebase.gov/catalog/item/5d949474e4b0c4f70d0db64a?format=json",
        "https://www.sciencebase.gov/catalog/item/5d949489e4b0c4f70d0db64d?format=json",
        "https://www.sciencebase.gov/catalog/item/5d94949de4b0c4f70d0db64f?format=json"
        ]

    # new dict
    files = {}

    for json_url in json_urls:
        try:
            r = requests.get(json_url, timeout=5)
            r.raise_for_status()
        except requests.exceptions.HTTPError as err_h:
            print("Http Error:", err_h)
            sys.exit(1)
        except requests.exceptions.ConnectionError as err_c:
            print("Error Connecting:", err_c)
            sys.exit(1)
        except requests.exceptions.Timeout as err_t:
            print("Timeout Error:", err_t)
            sys.exit(1)
        except requests.exceptions.RequestException as err:
            print(err)
            sys.exit(1)

        # load json
        tmp_json = json.loads(r.text)

        # loop & extract ancillary zip files
        for f in tmp_json['files']:
            if f['contentType'] == 'application/zip':
                f_dir = project_dir / 'ancillary_data'
                file = f_dir / f['name']
                shp = f_dir / f'{file.stem}.shp'
                # update the new dict
                files[f['title']] = {
                    'fname': f['name'],
                    'download': f['downloadUri'],
                    'f_dir': f_dir,
                    'file': file,
                    'shp': shp
                    }

    return files


def get_data(args):
    """
    Pass url to download the file through stream and unzip
    """
    v = args
    url, f_dir, file, shp = v['download'], v['f_dir'], v['file'], v['shp']

    n = file.name

    # create folders:
    if not f_dir.is_dir():
        os.makedirs(f_dir)

    if shp.is_file():
        # print(shp, ': exists!')
        return

    try:
        with requests.get(url, stream=True) as r:
            with open(file, 'wb') as f:
                shutil.copyfileobj(r.raw, f)
        print('Download Successful: ', n)
    except:
        print(f'Failed: {n}. Unable to retrieve the file!')

    try:
        archive = zipfile.ZipFile(file)
        archive.extractall(f_dir)
        print('Extract Successful: ', n)
    except:
        print('Failed to extract: ', n)


def validate_TauDEM():
    """
    run TauDEM with incomplete/missing args to retrieve
    an error message to validate install. Replace with more
    substantial test later
    """
    try:
        msg_bytes = subprocess.check_output('mpiexec -n 2 AreaD8')
        #bytes to string + slice 'Error'
        msg_str = msg_bytes.decode()[0:5]
        if msg_str == 'Error':
            print('Call to TauDEM via mpiexec successfully')
            return True
    except:
        print('Call to TauDEM via mpiexec failed. Please recheck and install TauDEM and its dependencies')
        sys.exit(1)


def update_config(ini, data_dict, project_dir):
    """update ini file with info from data dict"""

    config_file = configparser.ConfigParser(
        comment_prefixes=';',
        inline_comment_prefixes=';',
        allow_no_value=True,
        delimiters=':'
        )

    config_file.read(Path(ini))

    # get CBW and DRB physio files
    CBW_Physio = list(project_dir.rglob('*physiographic_regions_CBW.shp'))[0]
    DRB_Physio = list(project_dir.rglob('*physiographic_regions_DRB.shp'))[0]

    # update paths and flags
    if validate_TauDEM() == True:
        config_file.set('pre process dem', 'taudem', str(True))
        config_file.set('pre process dem', 'taudem cores', str(mp.cpu_count()))
    else:
        config_file.set('pre process dem', 'taudem', str(False))
    if CBW_Physio.is_file():
        config_file.set('file paths', 'physio cbw', str(CBW_Physio))
    if DRB_Physio.is_file():
        config_file.set('file paths', 'physio drb', str(DRB_Physio))
    if data_dict['census_roads']['shp'].is_file():
        config_file.set('file paths', 'census roads', str(data_dict['census_roads']['shp']))
    if data_dict['census_rails']['shp'].is_file():
        config_file.set('file paths', 'census rails', str(data_dict['census_rails']['shp']))

    # update ancillary dir
    ancil_dir = list(project_dir.rglob('*ancillary_data'))[0]
    if ancil_dir.is_dir():
        config_file.set('file paths', 'ancillary dir', str(ancil_dir))

    with open(ini, 'w') as cfile:
        config_file.write(cfile)

    return config_file


def main():

    intro_msg = '''
    +-+-+-+-+-+ +-+-+-+-+-+
    |F|A|C|E|T| |S|E|T|U|P|
    +-+-+-+-+-+ +-+-+-+-+-+

    setup.py script will automatically download ancillary datasets needed by FACET from Science Base.
    It also updates the config file with ancillary data file paths as well.
    The script takes anywhere from 30 minutes to 1-2 hours to run successfully depending on your internet speed.\n
    '''
    start = timer()
    # arg parser
    parser = argparse.ArgumentParser(
        prog='FACET set-up script',
        description=intro_msg
    )
    # user arg defined project path
    parser.add_argument('--path', '-p', help='full path for project', required=True)
    args = parser.parse_args()

    # intro message
    print(intro_msg)

    # project pathname
    project_dir = Path(args.path)

    # config file path
    config_ini_file = config.get_config_path()

    # get science base json data
    sb_json = json_to_dict(project_dir)

    ### STEP 1
    # download, extract and clean up
    st_step_01 = timer()
    print('''
    Step #1: Downloading and extracting ancillary datasets. This will take 30-45 mins. Please wait.
    ''')
    # number of items to download at any given time
    no_of_items = int(len(sb_json)/3)
    pool = mp.Pool(processes=no_of_items)
    pool.map(get_data, [(v) for v in sb_json.values()])
    pool.close()
    end_step_01 = round((timer()-st_step_01)/60.0, 2)
    print(f'Step #1 complete. Time elapsed: {end_step_01} mins\n')

    ### STEP 2
    st_step_02 = timer()
    # update and save config file
    print('Step #2: Updating config file and deleting temp files.')
    ini_dict = update_config(config_ini_file, sb_json, project_dir)

    # clean up zip files
    for v in sb_json.values():
        if v['file'].is_file():
            os.remove(v['file'])

    # result output for users
    end_step_02 = round((timer()-st_step_02)/60.0, 2)
    print(f'Step #2 complete. Time elapsed: {end_step_02} mins\n')

    files_downloaded = list(project_dir.rglob('*/ancillary_data/*.shp'))

    end = round((timer()-start)/60.0, 2)

    # message 1
    msg_1 = '''

    Setup complete. See summary below:

    1. Files downloaded:
    '''
    # message 2
    msg_2 = f'''
    2. Following variables were updated in config file:
        taudem: {ini_dict['pre process dem']['taudem']}
        taudem: {ini_dict['pre process dem']['taudem cores']}
        physio cbw: {ini_dict['file paths']['physio cbw']}
        physio drb: {ini_dict['file paths']['physio drb']}
        census roads: {ini_dict['file paths']['census roads']}
        census rails: {ini_dict['file paths']['census rails']}

        *************************************************
        *** CAUTION: READ THE SECTION BELOW CAREFULLY ***
        *************************************************

        Users need to open config.ini, review and edit the following parameters.
        For more details please refer to README, it will provide a detailed overview of each
        parameter and if it can be modified or not:

        * data_dir: This is a sub-directory inside your project directory where your data is stored
        * spatial ref: needs to be in PROJ.4 string (see README)
        * Under 'pre process dem':
            - taudem cores (default is set to using all cores, but users can manually change it)
            - resample resolution (see README)
        * breach options: 'pre-condition dem & fast-breach' is the default option (see README)
        * post processing: 'r exe path' needs to be changed based on your local installation of R

        DO NOT update values under the following subsections in config.ini (see README):
        * reach and order
        * cross section method
        * width from curvature via buff. method

    Total time elapsed: {end} mins'''

    for i in [msg_1, files_downloaded, msg_2]:
        if isinstance(i, str):
            print(i)
        if isinstance(i, list):
            for x in i:
                print(f'\t{x}')

if __name__ == '__main__':
    main()
