import pickle
import pandas as pd

def read_psmlist(path):
    with open(path, 'rb') as f:
        psmlist = pickle.load(f)
    return psmlist
    
def read_features(path):
    return pd.read_parquet(path)
