import logging
from pathlib import Path
import subprocess
import rasterio
import geopandas as gpd
import fsspec


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


def initialize_logger(log_file: str) -> logging.getLogger():
    """
    Initialize logger values and get logger object
    Args:
        log_file: Path to the log file

    Returns: Logger instance
    """

    clear_out_logger()

    if log_file.exists():
        log_file.unlink()

    logger = logging.getLogger("logger_loader")
    logging.basicConfig(filename=log_file, filemode="a")
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s [%(lineno)d] - %(message)s", "%m/%d/%Y %I:%M:%S %p"
    )
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


def clear_out_logger():
    """Remove all handlers associated with the root logger object"""
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)


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
