import geopandas as gpd
import numpy as np
from utils import utils

def create_flowline_qc_mask(flowlines, buffer, watershed):

    flowline_buffer = utils.vector_to_geodataframe(flowlines)
    watershed = utils.vector_to_geodataframe(watershed)
    flowline_buffer['geometry'] = flowline_buffer.buffer(buffer)
    flowline_buffer['id'] = 1
    flowline_buffer = flowline_buffer.dissolve(by='id')

    flowline_buffer = gpd.overlay(watershed, flowline_buffer, how='symmetric_difference')
    # flowline_buffer = flowline_buffer.drop(columns=['fdate'])
    return flowline_buffer


def create_waterbody_qc_mask(waterbody, ftypes, watershed):
    # [390, 436]
    waterbody = utils.vector_to_geodataframe(waterbody)
    watershed = utils.vector_to_geodataframe(watershed)
    waterbody = gpd.clip(waterbody, watershed)
    waterbody = waterbody.loc[waterbody['ftype'].isin(ftypes)]
    waterbody = gpd.overlay(waterbody, watershed, how='intersection')
    return waterbody



def flag_features_by_qc_mask(features, mask, flag, by_field="", output=""):
    uid = "uid"
    features[uid] = np.arange(features.shape[0])

    features_in_mask = gpd.sjoin(features, mask, how='inner', predicate='intersects')

    if by_field:
        uniques = set(features_in_mask[by_field])
        features.loc[features[by_field].isin(uniques), [flag]] = int(1)
    else:
        uniques = set(features_in_mask[uid])
        features.loc[features[uid].isin(uniques), [flag]] = int(1)

    if output:
        features.to_file(output)
    else:
        return features


