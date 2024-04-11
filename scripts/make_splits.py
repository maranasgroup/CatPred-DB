import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools, Draw
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import rdMolDescriptors as rdmd
from rdkit.Chem import Descriptors, rdChemReactions
import seaborn as sns
from tqdm import tqdm
import os
from sklearn.model_selection import GroupKFold, train_test_split

def assign_seq_clusters(df, id_cutoff):
    seqcol = 'sequence'
    import tempfile
    import uuid

    # Create a temporary directory and get a unique filename within that directory
    # with tempfile.TemporaryDirectory() as unique_filedir:
    unique_filedir = str(uuid.uuid4())
    # Using UUID for generating a unique filename
    unique_filename = unique_filedir + '_temp_seqs.fasta'

    print('Clustering seqs..')
    seqs_written = set()
    seqind_to_seq = {}
    f = open(unique_filename,'w')
    for ind,seq in enumerate(df[seqcol]):
        if pd.isna(seq): continue
        if not seq in seqs_written: 
            f.write(f'>{ind}\n{seq}\n')
            seqs_written.add(seq)
            seqind_to_seq[ind] = seq
        else: continue
    f.close()

    cmd = f'mmseqs easy-cluster {unique_filename} {unique_filedir}_clusterRes /tmp/ --min-seq-id {id_cutoff} -v 0'
    status = os.system(cmd)
    print(status, cmd)
    cluster_df = pd.read_csv(f'{unique_filedir}_clusterRes_cluster.tsv',sep='\t',names=['query_ind','target_ind'])
    os.system(f'rm -r {unique_filedir}*')
        
    q_to_tinds = {}
    for qind, tind in zip(cluster_df.query_ind,cluster_df.target_ind):
        if not qind in q_to_tinds:
            q_to_tinds[qind] = [tind]
        else: 
            q_to_tinds[qind].append(tind)
    c = 0
    ind_to_cluster = {}
    for q, ts in q_to_tinds.items():
        ind_to_cluster[q] = c
        for t in ts: ind_to_cluster[t] = c
        c+=1
    seq_to_cluster = {seqind_to_seq[ind]:c for ind, c in ind_to_cluster.items()}
    
    cluster_col = []
    for seq in df[seqcol]:
        if seq in seq_to_cluster:
            cluster_col.append(seq_to_cluster[seq])
        else:
            # print('no cluster found for ', seq)
            cluster_col.append(-1)
    df[f'{seqcol}_{int(100*id_cutoff)}cluster'] = cluster_col
    return df

import ipdb
def make_splits(df, frac = 0.1):    
    # ipdb.set_trace()
    df.reset_index(drop=True, inplace=True)
    train_df, test_df = train_test_split(df, test_size = frac, random_state=0)
    train_df.reset_index(drop=True, inplace=True)
    train_df, val_df = train_test_split(train_df, test_size = frac, random_state=0)
    train_df.reset_index(drop=True, inplace=True)
    val_df.reset_index(drop=True, inplace=True)
    test_df.reset_index(drop=True, inplace=True)
    
    print(f'Train: {len(train_df)*100/len(df)}')
    print(f'Val: {len(val_df)*100/len(df)}')
    print(f'Test: {len(test_df)*100/len(df)}')
    return train_df, val_df, test_df

def print_stats(train, test, parameter, dire, sim_list = [0.4, 0.6, 0.8, 0.99]):
    sequence_column = 'sequence'
    smiles_column = 'smiles'
    if parameter=='kcat': smiles_column='reactant_'+smiles_column
    else: smiles_column='substrate_'+smiles_column
    sim_list.reverse()
    print('-'*20)
    print('**Sequence stats**')
    print('-'*20)
    for sim in sim_list:
        simperc = int(sim*100)
        colname = sequence_column+f'_{simperc}cluster'
        print('At sim-cutoff:', sim,)
        now = test[~test[colname].isin(train[colname])]
        now.to_csv(f'{dire}/{parameter}-seq_test_{colname}.csv')
        perc = len(now)/len(test)
        print('Test perc.:', perc)
        print('Test number:', int(perc*len(test)))
        print('\n')
    print('-'*20)
        
def assign_all_clusters(df, is_kcat, sim_list = [0.4, 0.6, 0.8, 0.99]):
    for sim in sim_list:
        print(f'Assigning clusters at sim-cutoff = {sim}')
        df = assign_seq_clusters(df, id_cutoff=sim)
    return df

def make_and_savesplits(df, param, dire):
    if param=='kcat': is_kcat = True
    else: is_kcat = False
    frac = 0.1
    train, val, test = make_splits(df, frac)
    seqpre = 'random'
    if not os.path.exists(dire):
        os.mkdir(dire)
    train.to_csv(dire+param+f'-{seqpre}_train.csv')
    val.to_csv(dire+param+f'-{seqpre}_val.csv')
    test.to_csv(dire+param+f'-{seqpre}_test.csv')
    trainval = pd.concat([train,val]).reset_index(drop=True)
    trainval.to_csv(dire+param+f'-{seqpre}_trainval.csv')
    trainvaltest = pd.concat([trainval,test]).reset_index(drop=True)
    trainvaltest.to_csv(dire+param+f'-{seqpre}_trainvaltest.csv')
    print_stats(trainval, test, param, dire)
    return train, val, test
    
def main(args):
    datafile, parameter, outdir = args.input_file, args.param, args.save_dir
    df = pd.read_csv(datafile)
    df = assign_all_clusters(df, is_kcat=(parameter=='kcat'))#, 
    df_train, df_val, df_test = make_and_savesplits(df, parameter, dire=outdir)
    
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--param", help="parameter name", 
                        type=str, required=True)
    parser.add_argument("--input_file", help="input_file", 
                        type=str, required=True)
    parser.add_argument("--save_dir", help="dir to save splits", 
                        type=str, required=True)

    args, unparsed = parser.parse_known_args()
    
    return args

args = parse_args()
main(args)
