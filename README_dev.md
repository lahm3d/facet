# Floodplain and Channel Evaluation Tool (FACET)

This is FACET repo for continuous development (CD) and the published FACET 0.1.0 code can be accessed [here](https://code.usgs.gov/water/facet). There will be large gaps between published and CD repositories in terms of environment, behavior and algorithms. Once we reach a satisfactory point, we will publish the code on `https://code.usgs.gov`.

The repo has two branches `main` and `dev`. The main branch is similar to published code to keep up with library/module and environment updates along with any deprecated behaviors, and the `dev` branch is for experimental features.

## Getting started with FACET 0.2

Install git or github desktop and clone 

        https://github.com/lahm3d/facet.git - b dev


Install miniconda3 or anaconda: https://docs.conda.io/projects/miniconda/en/latest/

        cd c:/folder-where-you-cloned-or-unzipped-facet-repo

There are two options to install taudem:

Try option #1 once or twice, if that fails then simply follow option #2. There seems to be gdal version mismatch issues between taudem and other geospatial libs.

Option #1: Go to `environment.yml` and un-comment `taudem=5.8.3`. Skip to install conda environment section.

Option #2: Download and install [TauDEM](http://hydrology.usu.edu/taudem/taudem5/downloads.html) version 5.3.7, including the TauDEM dependencies.

Install conda environment

        conda env create -f environment.yml


### Download sample data: `https://gis-data.chesapeakebay.net/facet_misc_data/draft.zip`

Navigate to `facet/src/config.toml` and edit the following values:

-`batch_csv` : custom csv file where you can enter huc-id numbers to process, and to skip. Use the `facet/batch.csv` template file and modify as needed

- `ancillary`: Input paths for all the ancillary data which have been converted to geoparquet and stored on s3. Use the default paths for the files.

- `data`: Location of your data folder where you will be organizing all your watershed folders. This folder needs to be created manually.
- `streams`: streams or flowlines if you are using existing stream network 
- `nhd_area`: NHD area
- `nhd_waterbody`: NHD waterbody
- `nhd_physiography`: Physiography file
- `census_roads`: Census 2023 roads
- `census_rails`: Census 2023 rails

### Sample data
The folder structure has been reorganized but not finalized. If no version is provided in the config then all the files under `ver-1` will be created in the main folder.
```
└── huc_020600031001
    ├── 020600031001_dem.tif
    ├── 020600031001_flowlines.shp
    ├── 020600031001_hs.tif
    ├── 020600031001_watershed.shp
    └── ver-1
        ├── 020600031001_rr_crossings_ver-01.shp
        ├── 020600031001_burn_crossings_ver-01.tif
        ├── 020600031001_denoise_ver-01.tif
        ├── 020600031001_initiation_pixels_ver-01.tif
        ├── 020600031001_breach_ver-01.tif
        ├── 020600031001_slope_grid_sd8_ver-01.tif
        ├── 020600031001_d8_fdir_point_ver-01.tif
        ├── 020600031001_area_grid_ad8_ip_ver-01.tif
        ├── 020600031001_area_grid_ad8_ver-01.tif
        ├── 020600031001_network_coords_ver-01.dat
        ├── 020600031001_network_ver-01.shp
        ├── 020600031001_network_tree_ver-01.dat
        ├── 020600031001_sub_watersheds_ver-01.tif
        ├── 020600031001_network_order_ver-01.tif
        ├── 020600031001_slope_grid_dinf_ver-01.tif
        ├── 020600031001_flow_dir_dinf_ver-01.tif
        ├── 020600031001_hand_ver-01.tif
        └── 020600031001_sub_watersheds_ver-01.shp
```

To run facet, first open miniconda window:

        cd c:/folder-where-you-cloned-or-unzipped-facet-repo
        conda activate facet
        python src/facet.py


## Reporting bugs
Please consider reporting bugs and asking questions on the code.usgs.gov/FACET Issues page: [https://code.usgs.gov/water/facet/issues](https://code.usgs.gov/water/facet/issues). Users can also email questions to facet@usgs.gov.

## Acknowledgements

Funding for FACET was provided by a grant from the William Penn Foundation Delaware Watershed Research Fund and the U.S. Geological Survey.

## References

<a name="1">1</a>: C. Gangodagamage et al., "Wavelet-Compressed Representation of Landscapes for Hydrologic and Geomorphologic Applications," in IEEE Geoscience and Remote Sensing Letters, vol. 13, no. 4, pp. 480-484, April 2016.
doi: 10.1109/LGRS.2015.2513011

## Disclaimer

This software has been approved for release by the U.S. Geological Survey (USGS). Although the software has been subjected to rigorous review, the USGS reserves the right to update the software as needed pursuant to further analysis and review. No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the software is released on condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from its authorized or unauthorized use.

## License

This software is licensed under [CC0 1.0](http://creativecommons.org/publicdomain/zero/1.0/) and is in the [public domain](https://en.wikipedia.org/wiki/Public_domain) because it contains materials that originally
came from the [U.S. Geological Survey (USGS)](https://www.usgs.gov/), an agency of the [United States Department of Interior](https://www.doi.gov/). For more.
information, see the [official USGS copyright policy](http://www.usgs.gov/visual-id/credit_usgs.html#copyright/).

![Creative Commons logo](http://i.creativecommons.org/p/zero/1.0/88x31.png)
