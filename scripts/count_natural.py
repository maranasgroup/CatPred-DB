import json
import pandas as pd
from tqdm import tqdm
from rdkit import Chem

natural_list = set(json.load(open('./all_natural_metabolite_names.json')))
metabolite_inchi_smiles_dic = pd.read_csv('../datasets/metabolite_inchi_smiles_brenda_pubchem.tsv',sep='\t',index_col='metabolite')
name_to_smiles = {}
for ind, row in tqdm(metabolite_inchi_smiles_dic.iterrows()):
    try:
        name_to_smiles[row.name] = Chem.MolToSmiles(Chem.MolFromSmiles(row.smiles))
    except:
        continue
    
natural_smi_list = [name_to_smiles[name] for name in natural_list if name in name_to_smiles]

for parameter in ['kcat','km','ki']:
    print(parameter)
    datafile = f'./data/processed/splits_wpdbs/{parameter}-random_trainvaltest.csv'

    if parameter=='kcat':
        subcol = 'reactant_smiles'
    else:
        subcol = 'substrate_smiles'
        
    df = pd.read_csv(datafile)
    count = 0
    for ind, row in tqdm(df.iterrows()):
        natural = True
        substrates = row[subcol].split('.')
        for sub in substrates:
            if not sub in natural_smi_list:
                natural = False
        if natural: count+=1
    
    print(count)
