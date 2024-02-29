import pandas as pd

def generate_processing_batch(batch_csv):    
    batch = pd.read_csv(batch_csv, dtype=str)
    return batch.process.to_list()
