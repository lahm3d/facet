
import tomli
from pathlib import Path
from dataclasses import dataclass, InitVar, field

def read_toml(file):
    """
    read_toml _summary_

    _extended_summary_

    Parameters
    ----------
    file : _type_
        _description_

    Returns
    -------
    _type_
        _description_

    Raises
    ------
    FileNotFoundError
        _description_
    """
    if file.is_file():
        with open(file, mode="rb") as fp:
            return tomli.load(fp)
    else:
        raise FileNotFoundError(f'{file}: File not found!')


@dataclass
class CreateConfig:
    config: dict = field(default_factory=dict)

    def __post_init__(self):
        for key, value in self.config.items():
            setattr(self, key, value)


@dataclass
class CreateFilepaths:
    folder: Path
    huc: str
    version: str
    paths: InitVar[dict]

    def __post_init__(self, paths):
        for key, value in paths.items():
            stem, suffix = value.split('.')
            if self.version == "":
                parent = Path(self.folder) / f"huc_{self.huc}"
                basename = f"{self.huc}_{stem}.{suffix}"
            else:
                parent = Path(self.folder) / f"huc_{self.huc}" / self.version
                basename = f"{self.huc}_{stem}_{self.version}.{suffix}"
            fpath = parent / basename         
            setattr(self, key, fpath)

    def __repr__(self):
        attr_list = [f"{attr}={getattr(self, attr)}" for attr in self.__dict__ if not attr.startswith('_')]
        return f"{self.__class__.__name__}({', '.join(attr_list)})"


def class_to_dict(object):
    return object.__dict__


def create_config_and_filepath_objects(config_toml, paths_toml, huc):
    config = read_toml(config_toml)
    CONFIG = CreateConfig(config=config)

    paths = read_toml(paths_toml)
    PATHS = CreateFilepaths(
        folder=CONFIG.ancillary['data'],
        huc=huc,
        version=CONFIG.debug['version'],
        paths=paths
    )

    return CONFIG, PATHS