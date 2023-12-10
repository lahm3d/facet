import pandas as pd
from pathlib import Path

def generate_processing_batch(batch_csv):
    hucs = []
    
    batch = pd.read_csv(batch_csv, dtype=str)
    process = batch.process.to_list()
    skip = batch.skip.dropna().to_list()

    for huc in process:
        if huc not in skip:
            hucs.append(huc)

    return hucs 
