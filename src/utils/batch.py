import pandas as pd
from pathlib import Path

def generate_processing_batch(batch_csv, skip_csv):
    hucs = []
    
    batch = pd.read_csv(batch_csv, header=None, dtype=str)
    batch = batch[0].to_list()

    skip = pd.read_csv(skip_csv, header=None, dtype=str)
    skip = skip[0].to_list()

    for huc in batch:
        if huc not in skip:
            hucs.append(huc)

    return hucs 


