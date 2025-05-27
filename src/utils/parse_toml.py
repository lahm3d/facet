
import tomllib
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
            return tomllib.load(fp)
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
    paths_mod: dict
    preprocess_flag: bool
    preprocess_version: str
    xsec_flag: bool
    xsec_folder: str
    xsec_version: str
    
    def __post_init__(self, paths):
        if ( ( self.preprocess_flag == True ) & ( self.xsec_flag == False ) ):
            inputs = [
                'flowlines', 'dem', 'watershed', 'hs', 'physiography'
                ]
            for key, value in paths.items():
                stem, suffix = value.split('.')
    
                if stem in inputs:
                    parent = Path(self.folder) / self.huc
                    basename = f"{self.huc}_{stem}.{suffix}"
                else:
                    if self.version == "":
                        parent = Path(self.folder) / self.huc
                        basename = f"{self.huc}_{stem}.{suffix}"
                    else:
                        parent = Path(self.folder) / self.huc / self.version
                        basename = f"{self.huc}_{stem}_{self.version}.{suffix}"
                fpath = parent / basename         
                setattr(self, key, fpath)
                
        elif ( ( self.preprocess_flag == False ) & ( self.xsec_flag == False ) ):
            inputs = [
                'flowlines', 'dem', 'watershed', 'hs', 'physiography'
                ]
            preprocess = list( self.paths_mod['pre-process'].keys() )
            for key, value in paths.items():
                stem, suffix = value.split('.')
    
                if stem in inputs:
                    parent = Path(self.folder) / self.huc
                    basename = f"{self.huc}_{stem}.{suffix}"
                elif key in preprocess:
                    parent = Path(self.folder) / self.huc / self.preprocess_version
                    basename = f"{self.huc}_{stem}_{self.preprocess_version}.{suffix}"
                else:
                    if self.version == "":
                        parent = Path(self.folder) / self.huc
                        basename = f"{self.huc}_{stem}.{suffix}"
                    else:
                        parent = Path(self.folder) / self.huc / self.version
                        basename = f"{self.huc}_{stem}_{self.version}.{suffix}"
                fpath = parent / basename         
                setattr(self, key, fpath)
        elif ( ( self.preprocess_flag == False ) & ( self.xsec_flag == True ) ):
            inputs = [
                'flowlines', 'dem', 'watershed', 'hs', 'physiography'
                ]
            preprocess = list( self.paths_mod['pre-process'].keys() )
            xsecs = ['channel_xns', 'floodplain_xns', 'xn_coordinates']
            for key, value in paths.items():
                stem, suffix = value.split('.')
    
                if stem in inputs:
                    parent = Path(self.folder) / self.huc
                    basename = f"{self.huc}_{stem}.{suffix}"
                elif key in preprocess:
                    parent = Path(self.folder) / self.huc / self.preprocess_version
                    basename = f"{self.huc}_{stem}_{self.preprocess_version}.{suffix}"
                elif key in xsecs:
                    parent = Path( self.xsec_folder ) / self.huc / self.xsec_version
                    basename = f"{self.huc}_{stem}_{self.xsec_version}.{suffix}"
                else:
                    if self.version == "":
                        parent = Path(self.folder) / self.huc
                        basename = f"{self.huc}_{stem}.{suffix}"
                    else:
                        parent = Path(self.folder) / self.huc / self.version
                        basename = f"{self.huc}_{stem}_{self.version}.{suffix}"
                fpath = parent / basename         
                setattr(self, key, fpath)            
        elif ( ( self.preprocess_flag == True ) & ( self.xsec_flag == True ) ):
            inputs = [
                'flowlines', 'dem', 'watershed', 'hs', 'physiography'
                ]
            xsecs = ['channel_xns', 'floodplain_xns', 'xn_coordinates']
            for key, value in paths.items():
                stem, suffix = value.split('.')
    
                if stem in inputs:
                    parent = Path(self.folder) / self.huc
                    basename = f"{self.huc}_{stem}.{suffix}"
                elif key in xsecs:
                    parent = Path( self.xsec_folder ) / self.huc / self.xsec_version
                    basename = f"{self.huc}_{stem}_{self.xsec_version}.{suffix}"
                else:
                    if self.version == "":
                        parent = Path(self.folder) / self.huc
                        basename = f"{self.huc}_{stem}.{suffix}"
                    else:
                        parent = Path(self.folder) / self.huc / self.version
                        basename = f"{self.huc}_{stem}_{self.version}.{suffix}"
                fpath = parent / basename         
                setattr(self, key, fpath)
        # add parent folder 
        setattr(self, "parent", Path(self.folder) / self.huc)


    def __repr__(self):
        attr_list = [f"{attr}={getattr(self, attr)}" for attr in self.__dict__ if not attr.startswith('_')]
        return f"{self.__class__.__name__}({', '.join(attr_list)})"


def class_to_dict(object):
    return object.__dict__


def create_config(config_toml):
    config = read_toml(config_toml)
    return CreateConfig(config=config)

def create_filepaths(paths_toml, config, huc):
    paths_mod = read_toml(paths_toml)
    paths = dict( paths_mod['inputs'].items() )
    paths.update( paths_mod['pre-process'].items() )
    paths.update( paths_mod['cross-section-method'].items() )
    Paths = CreateFilepaths(
        folder = config.ancillary['data'],
        huc = huc,
        version = config.debug['version'],
        paths = paths,
        paths_mod = paths_mod, 
        preprocess_flag = config.preprocess['preprocess_flag'],
        preprocess_version = config.preprocess['preprocess_version'],
        xsec_flag = config.reuse_xn['reuse_xn_flag'],
        xsec_folder = config.reuse_xn['reuse_xn_data'],
        xsec_version = config.reuse_xn['reuse_xn_version']
    )

    return Paths
