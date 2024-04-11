import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
import ipdb
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--brenda_file", help="full path to brenda processed data file", type=str, required=True)
    parser.add_argument("--sabio_file", help="full path to sabio processed data file", type=str, required=True)
    parser.add_argument("--output_dir", help="full path to output directory to save merged data files", type=str, required=True)
    args, unparsed = parser.parse_known_args()
    parser = argparse.ArgumentParser()

    return args

args = parse_args()

sabio_df = pd.read_csv(f'{args.sabio_file}')
brenda_df = pd.read_csv(f'{args.brenda_file}')

# reaction_smiles can be not sorted and appearing as non unique
# For example A.B smiles is not same as B.A
# By sorting they can be made equal
reacsmi_col = []
for ind, row in sabio_df.iterrows():
    smi = row.reaction_smiles
    if not pd.isna(smi): 
        reacs, prods = smi.split('>>')
        reacs = sorted(reacs.split('.'))
        prods = sorted(prods.split('.'))
        smi = ".".join(reacs) + '>>' + ".".join(prods)
    reacsmi_col.append(smi)
sabio_df['reaction_smiles'] = reacsmi_col

reacsmi_col = []
for ind, row in brenda_df.iterrows():
    smi = row.reaction_smiles
    if not pd.isna(smi): 
        reacs, prods = smi.split('>>')
        reacs = sorted(reacs.split('.'))
        prods = sorted(prods.split('.'))
        smi = ".".join(reacs) + '>>' + ".".join(prods)
    reacsmi_col.append(smi)
brenda_df['reaction_smiles'] = reacsmi_col

cols = ['sequence','sequence_source','uniprot','reaction_smiles','value','temperature','ph','taxonomy_id']#

kcat_brenda = brenda_df[brenda_df.parameter=='turnover_number']
kcat_brenda = kcat_brenda[cols]
kcat_brenda.reset_index(inplace=True,drop=True)
    
kcat_sabio = sabio_df[sabio_df.param_name=='kcat']
kcat_sabio['value'] = kcat_sabio.param_value_stdunit
kcat_sabio = kcat_sabio[cols]

kcat_merged = pd.concat([kcat_brenda,kcat_sabio])
values = []
for val in kcat_merged.value:
    if val>1e6: val=1e6
    elif val<1e-6: val=1e-6
    values.append(val)
kcat_merged['value'] = values
# cap very large or very small values

kcat_merged['log10_value'] = np.log10(kcat_merged.value)
kcat_merged.drop_duplicates(inplace=True) 
kcat_merged[['reactant_smiles','product_smiles']] = kcat_merged['reaction_smiles'].str.split('>>',expand=True)
kcat_merged.dropna(subset=['ec','taxonomy_id','sequence','reactant_smiles','log10_value'],inplace=True)
kcat_merged.reset_index(inplace=True,drop=True)

# km
cols = ['sequence','sequence_source','uniprot','substrate_smiles','value','temperature','ph','taxonomy_id']

km_brenda = brenda_df[brenda_df.parameter=='km_value']
km_brenda = km_brenda[cols]

km_sabio = sabio_df[sabio_df.param_name=='Km']
km_sabio['value'] = km_sabio.param_value_stdunit
km_sabio['substrate_smiles'] = km_sabio.param_species_smiles
km_sabio = km_sabio[cols]

km_merged = pd.concat([km_brenda,km_sabio])
values = []
for val in km_merged.value:
    if val>1e4: val=1e4
    elif val<1e-8: val=1e-8
    values.append(val)
km_merged['value'] = values

km_merged['log10_value'] = np.log10(km_merged.value)
km_merged.drop_duplicates(inplace=True) 
km_merged.dropna(subset=['ec','taxonomy_id','sequence','substrate_smiles','log10_value'],inplace=True)
km_merged.reset_index(inplace=True,drop=True)

# ki
cols = ['sequence','sequence_source','uniprot','substrate_smiles','value','temperature','ph','taxonomy_id']

ki_brenda = brenda_df[brenda_df.parameter=='ki_value']
ki_brenda = ki_brenda[cols]

ki_sabio = sabio_df[sabio_df.param_name=='Ki']
ki_sabio['value'] = ki_sabio.param_value_stdunit
ki_sabio['substrate_smiles'] = ki_sabio.param_species_smiles
ki_sabio = ki_sabio[cols]

ki_merged = pd.concat([ki_brenda,ki_sabio])
values = []
for val in ki_merged.value:
    if val>1e4: val=1e4
    elif val<1e-10: val=1e-10
    values.append(val)
ki_merged['value'] = values

ki_merged['log10_value'] = np.log10(ki_merged.value)
ki_merged.drop_duplicates(inplace=True) 
ki_merged.dropna(subset=['ec','taxonomy_id','sequence','substrate_smiles','log10_value'],inplace=True)
ki_merged.reset_index(inplace=True,drop=True)

def handle_duplicates(df, param):
    if param=='kcat': grouper = ['reactant_smiles','sequence']
    else: grouper = ['substrate_smiles','sequence']
    newdf = pd.DataFrame()
    groups = df.groupby(grouper)
    for each in tqdm(groups):
        groupname, group = each
        if param=='kcat':
            value = group.log10_value.max()
            group['log10kcat_max'] = value
        else: 
            value = group.log10_value.mean()
            group[f'log10{param}_mean'] = value
        
        group['group'] = '__'.join(groupname)
        newdf = pd.concat([newdf,group.iloc[:1]])
    return newdf.reset_index(drop=True)

ipdb.set_trace()

kcat_merged = handle_duplicates(kcat_merged,'kcat')
km_merged = handle_duplicates(km_merged,'km')
ki_merged = handle_duplicates(ki_merged,'ki')

if not os.path.exists(args.output_dir):
    os.mkdir(args.output_dir)

kcat_merged.to_csv(f'{args.output_dir}/kcat_merged_max.csv')
km_merged.to_csv(f'{args.output_dir}/km_merged_mean.csv')
ki_merged.to_csv(f'{args.output_dir}/ki_merged_mean.csv')

print(len(kcat_merged))
print(len(km_merged))
print(len(ki_merged))

import sys
sys.exit(0)

from sklearn.model_selection import train_test_split

def split_train_test(df):
    train_df, test_df = train_test_split(df, test_size=0.1, random_state=0)
    train_df, val_df = train_test_split(train_df, test_size=0.1, random_state=0)
    train_df['split'] = 'train'
    val_df['split'] = 'val'
    test_df['split'] = 'test'
    
    return pd.concat([train_df, val_df]).reset_index(drop=True), test_df


kcat_train, kcat_test = split_train_test(kcat_merged)
km_train, km_test = split_train_test(km_merged)
ki_train, ki_test = split_train_test(ki_merged)

kcat_train.to_csv('./final_data/kcat_train.csv')
kcat_test.to_csv('./final_data/kcat_test.csv')
km_train.to_csv('./final_data/km_train.csv')
km_test.to_csv('./final_data/km_test.csv')
ki_train.to_csv('./final_data/ki_train.csv')
ki_test.to_csv('./final_data/ki_test.csv')

kcat_train_seq = kcat_train.dropna(subset=['sequence']).reset_index(drop=True)
kcat_test_seq = kcat_test.dropna(subset=['sequence']).reset_index(drop=True)
km_train_seq = km_train.dropna(subset=['sequence']).reset_index(drop=True)
km_test_seq = km_test.dropna(subset=['sequence']).reset_index(drop=True)
ki_train_seq = ki_train.dropna(subset=['sequence']).reset_index(drop=True)
ki_test_seq = ki_test.dropna(subset=['sequence']).reset_index(drop=True)

# stats for data visualization

# 
ecs = set(kcat_merged.ec).union(set(km_merged.ec)).union(set(ki_merged.ec))
ecs_seq = set(kcat_train_seq.ec).union(set(kcat_test_seq.ec)).union(set(km_train_seq.ec)).union(set(km_test_seq.ec)).union(set(ki_train_seq.ec)).union(set(ki_test_seq.ec))

ec_to_seq = {ec:0 for ec in ecs}
for ec in ecs_seq: ec_to_seq[ec] = 1

smis = set(km_merged.ec).union(set(ki_merged.ec))
for reacsmi in kcat_merged.reactant_smiles:
    reacs = smis.split('.')
    for r in reacs:smis.add(r)

smis_seq = set(km_train_seq.ec).union(set(km_test_seq.ec)).union(set(ki_train_seq.ec)).union(set(ki_test_seq.ec))
for reacsmi in kcat_train_seq.reactant_smiles:
    reacs = smis_seq.split('.')
    for r in reacs:smis_seq.add(r)
for reacsmi in kcat_test_seq.reactant_smiles:
    reacs = smis.split('.')
    for r in reacs:smis_seq.add(r)
    
smi_to_seq = {smi:0 for smi in smis}
for smi in smis_seq: smi_to_seq[smi] = 1


