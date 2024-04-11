import requests
from joblib import Parallel, delayed
import os
import argparse
import sys
import json
import pandas as pd
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_file", help="data file", type=str, required=True)
    parser.add_argument("--alphafold_dir", help="directory where alphafold models are stroed", type=str, required=True)
    parser.add_argument("--esm_dir", help="directory where to store esm models", type=str, required=True)
    parser.add_argument("--out_file", help="out file", type=str, required=True)
    
    args, unparsed = parser.parse_known_args()
    parser = argparse.ArgumentParser()

    return args

args = parse_args()

def _is_empty_model(modelpath):
    try:
        structure = esm.inverse_folding.util.load_structure(modelpath, 'A')
        return False
    except:
        return True

import ipdb
ipdb.set_trace()
alphafold_models = os.listdir(args.alphafold_dir)
esm_models = os.listdir(args.esm_dir)

af_model_path = lambda uni: args.alphafold_dir + f'/AF-{uni}-F1-model_v4.pdb'
esm_model_path = lambda uni: args.esm_dir + f'/ESMFold-{uni}-v1.pdb'

df = pd.read_csv(args.data_file)

modelpaths = []
for uni in tqdm(df.uniprot):
    modelpath = af_model_path(uni)
    modelpath2 = esm_model_path(uni)
    if not _is_empty_model(modelpath): modelpaths.append(modelpath)
    elif not _is_empty_model(modelpath2): modelpaths.append(modelpath2)
    else: modelpaths.append(None)
df['pdbpath'] = modelpaths

df.info()
df.dropna(subset=['pdbpath'],inplace=True)
df.reset_index(inplace=True,drop=True)
df.to_csv(args.out_file)