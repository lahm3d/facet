import logging
import subprocess
from pathlib import Path
from time import perf_counter

import geopandas as gpd
import rasterio


def elapsed_time(start):
    elapsed = perf_counter() - start  # elapsed time in seconds (float)

    hours = elapsed / 3600
    minutes = elapsed / 60
    seconds = elapsed

    if hours >= 1:
        return f'{hours:.2f} hrs'
    elif minutes >= 1:
        return f'{minutes:.2f} mins'
    else:
        return f'{seconds:.1f} seconds'

def create_folder(Paths):

    if Paths.version != "":
        (Path(Paths.parent) / Paths.version).mkdir(parents=True, exist_ok=True)
    else:
        (Path(Paths.parent) / Paths.version).mkdir(parents=True, exist_ok=True)


def run_command(cmd: str, logger: logging.getLogger()) -> None:
    """
    Execute commands as subprocesses

    Args:
        cmd: Command to run as a string
        logger: Logger instance

    Returns: None
    """
    try:
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        output, err = p.communicate()

        # Get some feedback from the process to print out:
        if err is None:
            text = output.decode()
            print("\n", text, "\n")
        else:
            print(err)

    except subprocess.CalledProcessError as e:
        logger.critical(f"failed to return code: {e}")
    except OSError as e:
        logger.critical(f"failed to execute shell: {e}")
    except IOError as e:
        logger.critical(f"failed to read file(s): {e}")


def vector_to_geodataframe(file):
    file = Path(file)
    ext = file.suffix

    if ext == '.parquet':
        with fsspec.open(file) as parquet:
            return gpd.read_parquet(parquet)
    else:
        return gpd.read_file(file)


class NoWarningsFilter(logging.Filter):
    def filter(self, record):
        # Return False to block WARNING messages, True to allow others
        return record.levelno != logging.WARNING


def initialize_logger(log_file: Path) -> logging.Logger:
    """
    Initialize logger values and get logger object
    Args:
        log_file: Path to the log file

    Returns: Logger instance
    """

    logger = logging.getLogger(__name__)

    # Clear existing handlers from this specific logger to prevent duplicates
    if logger.hasHandlers():
        logger.handlers.clear()

    # log formatting
    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s() [%(lineno)d] --> %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    # Add FileHandler
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Stream handler logs INFO and ERROR but excludes WARNING via filter
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)
    stream_handler.addFilter(NoWarningsFilter())
    logger.addHandler(stream_handler)

    logger.setLevel(logging.DEBUG)

    return logger


def my_callback(value):
    if not "%" in value:
        print(value)


def rasterize(vector_file: Path, ref_raster: Path, raster_file: Path, ID: str) -> None:
    """
    Rasterizes a geodataframe based on a template raster
    Args:
        vector_file: Path to the vector file to rasterize
        ref_raster: Path to the template raster
        raster_file: Path to the output raster

    Returns: None, writes the rasterized file to disk

    Thanks to:
    https://gis.stackexchange.com/questions/151339/rasterize-a-shapefile-with-geopandas-or-fiona-python
    """
    gdf = gpd.read_file(vector_file)

    with rasterio.open(ref_raster) as rst:
        meta = rst.meta.copy()
    meta.update(compress="lzw")
    meta.update(dtype=rasterio.int32)
    meta.update(nodata=0)

    with rasterio.open(raster_file, "w+", **meta) as out:
        out_arr = out.read(1)
        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom, value) for geom, value in zip(gdf.geometry, gdf[ID]))
        arr_burned = rasterio.features.rasterize(
            shapes=shapes, fill=0, out=out_arr, transform=out.transform
        )
        out.write_band(1, arr_burned)