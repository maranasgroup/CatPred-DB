import os
import json
import time
import pandas as pd
import numpy as np
import requests
import ete3
import re
import rdkit
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import rdChemReactions
from joblib import delayed, Parallel
from rdkit.Chem import Descriptors
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import requests
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rawdata_dir", help="raw data directory", type=str, required=True)
    parser.add_argument("--processed_dir", help="processed dir to save", 
                        type=str, required=True)
    parser.add_argument("--out_file_name", help="name for out file", 
                        type=str, required=True)

    args, unparsed = parser.parse_known_args()
    parser = argparse.ArgumentParser()

    return args

args = parse_args()
PARAMETERS_OF_INTEREST = ['km_value','kcat_km','ki_value','ic50','turnover_number']
EC_WORD_NAMES = ["ec1", "ec2", "ec3", "ec"]
TAX_WORD_NAMES = [
"superkingdom",
"phylum",
"class",
"order",
"family",
"genus",
"species"]
INCLUDE_OTHERS_AS_WILDTYPE = True
N_PROCS = 30

DATA_DIR = args.rawdata_dir
out_path = args.processed_dir + '/'+ args.out_file_name

# Load BRENDA raw data
# Downloaded from https://www.brenda-enzymes.org/download.php
data = json.load(open(f'{DATA_DIR}/brenda_2022_2.json'))['data']

## Ligand data
# Go to https://www.brenda-enzymes.org/search_result.php?a=13 and place a blank query to get all ligands and their InChi strings
# Downloaded and pre-processed this along with missing metabolites using PubChem id-exchange service
# https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi

metabolite_inchi_smiles_dic = pd.read_csv(f'{DATA_DIR}/metabolite_inchi_smiles_brenda_pubchem.tsv',sep='\t',
                                          index_col='metabolite')

# Create an empty dataframe and columns
df = pd.DataFrame()

# storing ec numbers
eccol = []
# storing organism names
orgcol = []

# storing parameter type (like turnover_number, km_value etc.)
paramcol = []
# storing parameter values (like km,ki etc.)
valcol = []

# storing reactions (reactants, products)
rxncol = []

# storing substrates
subcol = []

# storing if the substrate is from a natural reaction or not
natural_substrate_col = []

# storing uniprots wherever available
unicol = []

# storing comment strings
commentcol = []

# storing ph optimum & temperature optimum
phopt_col = []
topt_col = []

# all
all_metabolite_names = set()
all_natural_metabolite_names = set()

#metals and ions list
metals_ions_list = set()

for ec in tqdm(data):
    if not ('proteins' in data[ec] or 'organisms' in data[ec]):
        # this means we cannot featurize entries from this ec, so skip it
        continue
    
    # not all entries for an ec number have all params
    # see which among the parameters of interest are present
    params_present = []
    for parameter_data in PARAMETERS_OF_INTEREST:
        if parameter_data in data[ec]: params_present.append(parameter_data)
    
    # for mapping organisms
    if 'organisms' in data[ec]: orgs_ec = data[ec]['organisms']
    else: orgs_ec = {}
    if 'proteins' in data[ec]: orgs_protein = data[ec]['proteins']
    else: orgs_protein = orgs_ec
    
    metals_ions_now = set()
    if 'metals_ions' in data[ec]:
        for metalion_dict in data[ec]['metals_ions']:
            metals_ions_now.add(metalion_dict['value'])
            metals_ions_list = metals_ions_list.union(metals_ions_now)
    
    cofactors_now = set()
    if 'cofactor' in data[ec]:
        for cof_dict in data[ec]['cofactor']:
            cofactors_now.add(cof_dict['value'])
    # for mapping optimum temperatures by organism
    topt_by_org = {}
    if 'temperature_optimum' in data[ec]:
        topt_data = data[ec]['temperature_optimum']
        for each in topt_data:
            if 'organisms' in each:
                org_topt = each['organisms'][0]
                if 'num_value' in each:
                    topt_val = each['num_value']
                elif 'min_value' in each and 'max_value' in each:
                    topt_val = (each['min_value'] + each['max_value'])/2.0
                else:
                    topt_val = None
                    continue
                    
                if org_topt in topt_by_org:
                    topt_by_org[org_topt].append(topt_val)
                else:
                    topt_by_org[org_topt] = [topt_val]
    
    # for mapping optimum ph by organism
    ph_by_org = {}
    if 'ph_optimum' in data[ec]:
        ph_data = data[ec]['ph_optimum']
        for each in ph_data:
            if 'organisms' in each:
                org_ph = each['organisms'][0]
                if 'num_value' in each:
                    ph_val = each['num_value']
                elif 'min_value' in each and 'max_value' in each:
                    ph_val = (each['min_value'] + each['max_value'])/2.0
                else:
                    ph_val = None
                    continue
                if org_ph in ph_by_org:
                    ph_by_org[org_ph].append(ph_val)
                else:
                    ph_by_org[org_ph] = [ph_val]
    

    # collect all the reactions of this EC grouped by substrates
    sub_to_reactions = {}
    
    def _reaction_string(reacs,prods):
        return f'{" + ".join(reacs)} >> {" + ".join(prods)}'
    
    if 'reaction' in data[ec]: 
        for rxn in data[ec]['reaction']:
            try:
                reacs = rxn['educts']
                prods = rxn['products']
            except KeyError:
                reacs = []
                prods = []
            
            each_rxn_str = _reaction_string(reacs,prods)
            
            for each in reacs:
                if each in sub_to_reactions:
                    sub_to_reactions[each].add(each_rxn_str)
                else:
                    sub_to_reactions[each] = set([each_rxn_str])
                
    # collect any reactant or product from natural_reaction (s) as natural_substrates
    natural_substrates = []
    if 'natural_reaction' in data[ec]: 
        for rxn in data[ec]['natural_reaction']:
            found = False
            try:
                reacs = rxn['educts']
                prods = rxn['products']
            except KeyError:
                reacs = []
                prods = []
                
            each_rxn_str = _reaction_string(reacs,prods)
            natural_substrates.extend(reacs)
            natural_substrates.extend(prods)
            for natural_met in reacs+prods:
                all_natural_metabolite_names.add(natural_met)
            for each in reacs:
                if each in sub_to_reactions:
                    sub_to_reactions[each].add(each_rxn_str)
                else:
                    sub_to_reactions[each] = set([each_rxn_str])
        
    # now collect parameters
    for parameter_data in params_present:
        for entry in data[ec][parameter_data]:
            if 'num_value' in entry and 'value' in entry and 'organisms' in entry:
                val = entry['num_value']
                if 'organisms' in entry:
                    org_now = entry['organisms']
                    prot_now = org_now
                elif 'proteins' in entry:
                    prot_now = entry['proteins']

                subname = entry['value']
                if parameter_data=='turnover_number' and (subname in cofactors_now or subname in metals_ions_now):
                    continue
                if not pd.isna(val) and not pd.isna(subname):
                    try:
                        val = float(val)
                    except ValueError:
                        continue
                    # get orgname if present, skip entry otherwise
                    try:
                        orgname = data[ec]['organisms'][org_now[0]]['value']
                    except KeyError:
                        continue
                    # get uniprot accessions , organism names

                    unis = []
                    orgnames = []

                    for org in prot_now:
                        if org in orgs_protein:
                            for each in orgs_protein[org]:
                                if 'accessions' in each:
                                    unis.append(each['accessions'])
                                else:
                                    unis.append(None)

                    for org in org_now:
                        if org in orgs_ec:
                            if 'value' in orgs_ec[org]: 
                                orgnames.append(orgs_ec[org]['value'])
                            else:
                                orgnames.append(None)

                    orgcol.append(orgnames[0])
                    if org_now[0] in ph_by_org:
                        phopt_col.append(ph_by_org[org_now[0]])
                    else:
                        phopt_col.append([])
                        
                    if org_now[0] in topt_by_org:
                        topt_col.append(topt_by_org[org_now[0]])
                    else:
                        topt_col.append([])
                        
                    unis_ = []
                    for u in unis:
                        if u is None: continue
                        unis_.append(u)

                    if not unis_: unis_ = None

                    unicol.append(unis_)
                    valcol.append(val)
                    paramcol.append(parameter_data)
                    subcol.append(subname)
                    all_metabolite_names.add(subname)
                    if subname in natural_substrates:
                        natural_substrate_col.append(True)
                    else:
                        natural_substrate_col.append(False)
                    eccol.append(ec)
                    rxncol.append(sub_to_reactions)
                    # natrxncol.append(natural_reactions)

                    if 'comment' in entry:
                        commentcol.append(entry['comment'])
                    else:
                        commentcol.append('')
                    
df['value'] = valcol
df['parameter'] = paramcol
df['substrate'] = subcol
df['natural_substrate'] = natural_substrate_col
df['organism'] = orgcol
df['ph_opt'] = phopt_col
df['temp_opt'] = topt_col
df['comment'] = commentcol
df['uniprot'] = unicol
df['ec'] = eccol
df['reactions'] = rxncol

import ipdb
ipdb.set_trace()

f = open('all_natural_metabolite_names.json','w')
f.write(json.dumps(list(all_natural_metabolite_names), indent=True))
f.close()

sys.exit(0)

# Get NCBI Taxonomy parser - to convert organism names into their taxonomic lineages
ncbi = ete3.NCBITaxa()
get_taxid_from_organism = lambda organism: ncbi.get_name_translator([organism])[organism][0]

ec_embed_cols = ["ec1", "ec2", "ec3", "ec"]
tax_embed_cols = [
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

def get_ec_words(ec):
    if '-' in ec: ec.replace('-','UNK')
    ec_chars = ec.split('.')
    ec_words = {f"ec{i}": '.'.join(ec_chars[:i]) for i in range(1,4)}
    ec_words['ec'] = ec
    return ec_words

def get_tax_words(taxid, ncbi):   
    try:
        lineage = ncbi.get_lineage(taxid)
        rank_dict = ncbi.get_rank(lineage)
        rank_dict_return = {}
        for rankid, rankname in rank_dict.items():
            if rankname.lower() in tax_embed_cols: rank_dict_return[rankname.lower()] = ncbi.get_taxid_translator([rankid])[rankid]
    except:
        rank_dict_return = {tax: 'UNK' for tax in tax_embed_cols}
    return rank_dict_return

taxid_col = []
for ind, row in tqdm(df.iterrows()):
    org = row.organism
    try:
        taxid = get_taxid_from_organism(org)
    except KeyError:
        taxid = None
    taxid_col.append(taxid)
df['taxonomy_id'] = taxid_col

ec_words = []
for ind, row in df.iterrows():
    words = get_ec_words(row.ec)
    ec_words.append(words)
for col in ec_embed_cols:
    col_values = [ec_words[i][col] for i in range(len(df))]
    df[col] = col_values

tax_words = []
for ind, row in df.iterrows():
    words = get_tax_words(row.taxonomy_id, ncbi)
    tax_words.append(words)
for col in tax_embed_cols:
    col_values = []
    for i in range(len(df)):
        if col in tax_words[i]:
            col_values.append(tax_words[i][col])
        else:
            col_values.append('UNK')
    df[col] = col_values

# import ipdb
# ipdb.set_trace()

# now SMILES mapping
smiles_col = []
for ind, row in tqdm(df.iterrows()):
    sub = row.substrate
    if not sub in metabolite_inchi_smiles_dic.index: 
        smiles_col.append(None)
    else: 
        smiles_col.append(metabolite_inchi_smiles_dic.loc[sub].smiles)
        
df['substrate_smiles'] = smiles_col

mwcol = []

for smi in tqdm(metabolite_inchi_smiles_dic.smiles):
    mw = None
    if not smi is None: 
        try:
            mol = Chem.MolFromSmiles(smi)
            mw = Descriptors.MolWt(mol)
        except:
            pass
    mwcol.append(mw)
metabolite_inchi_smiles_dic['MW'] = mwcol

def sort_by_second(item): return item[1]

smi_to_mw_dic = {}
for smi, mw in zip(metabolite_inchi_smiles_dic.smiles, metabolite_inchi_smiles_dic.MW):
    smi_to_mw_dic[smi] = mw
    
def sum_mwts(smis):
    mw = 0
    for smi in smis:
        mw+=smi_to_mw_dic[smi]
    return mw

# Assign enzyme type by parsing comments
# Assign everything as wild type first and filter out entries that are 'not possibly' wild-type 
df["enzyme_type"] = np.nan
df.loc[pd.isnull(df["comment"])] = ""
df["enzyme_type"] = "wild type"
df["enzyme_type"][df['comment'].str.contains("mutant")] = "mutant"
df["enzyme_type"][df['comment'].str.contains("mutate")] = "mutant"
df["enzyme_type"][df['comment'].str.contains("chimera")] = "mutant"
df["enzyme_type"][df['comment'].str.contains("inhibitor")] = "inhibition"
df["enzyme_type"][df['comment'].str.contains("inhibition")] = "inhibition"
df["enzyme_type"][df['comment'].str.contains("presence of")] = "regulated"
df["enzyme_type"][df['comment'].str.contains("recombinant")] = "recombinant"
df["enzyme_type"][df['comment'].str.contains("allozyme")] = "allozyme"
df["enzyme_type"][df['comment'].str.contains("alloenzyme")] = "allozyme"
df["enzyme_type"][df['comment'].str.contains("isozyme")] = "isozyme"
df["enzyme_type"][df['comment'].str.contains("isoenzyme")] = "isozyme"
df["enzyme_type"][df['comment'].str.contains("isoform")] = "isozyme"

# some entries still belong to mutants 
# check for X__Y type of comments (eg: W358A)
# these correspond to mutations should be removed
pat = r'[ACDEFGHIKLMNPQRSTVWY][0-9]+[ACDEFGHIKLMNPQRSTVWY]'
mutations_col = []
enztype_col = []
for com,enz_type in zip(df.comment,df.enzyme_type):
    items = re.findall(pattern=pat, string=com)
    if len(items)>0:
        mutations_col.append(';'.join(items))
        enztype_col.append('mutant')
    else:
        mutations_col.append(None)
        enztype_col.append(enz_type)
df['enzyme_type'] = enztype_col
df['mutations'] = mutations_col

# Now assigne ph and temperature first from comments
# then from optimum values for organism, ec
phrow = []
temprow = []
for com,phopt,topt in zip(df.comment,df.ph_opt,df.temp_opt):
    pat = r'[0-9][0-9].C'
    try:
        temp = int(re.findall(pattern=pat, string=com)[0].split('C')[0][:2])
        pat = r'pH [0-9]\.[0-9]'
        ph = float(re.findall(pattern=pat, string=com)[0].split('pH')[-1].strip())
    except:
        if phopt: ph = np.average(phopt)
        else: ph = None
        if topt: temp = np.average(topt)
        else: temp = None
    phrow.append(ph)
    temprow.append(temp)
    
df['ph'] = phrow
df['temperature'] = temprow

# fill in for missing ph, temperature using org,ec groups
orgec_to_topt = {}
orgec_to_phopt = {}
org_to_topt = {}
org_to_phopt = {}
for _, row in df.iterrows():
    org = row.organism
    ec = row.ec
    orgec = org+'__'+ec
    if not orgec in orgec_to_topt:
        orgec_to_topt[orgec] = []
    if not orgec in orgec_to_phopt:
        orgec_to_phopt[orgec] = []
    if not org in org_to_topt:
        org_to_topt[org] = []
    if not org in org_to_phopt:
        org_to_phopt[org] = []
    if not pd.isna(row.temperature): 
        orgec_to_topt[orgec].append(row.temperature)
        org_to_topt[org].append(row.temperature)
    if not pd.isna(row.ph): 
        orgec_to_phopt[orgec].append(row.ph)
        org_to_phopt[org].append(row.ph)

tempcol = []
phcol = []
for org, ec, temp, ph in zip(df.organism, df.ec, df.temperature,df.ph):
    orgec = org+'__'+ec
    if pd.isna(temp):
        if len(orgec_to_topt[orgec])>0:
            temp = np.median(orgec_to_topt[orgec])
        else:
            if len(org_to_topt[org])>0:
                temp = np.median(org_to_topt[org])
            else:
                temp = None
    if pd.isna(ph):
        if len(orgec_to_phopt[orgec])>0:
            ph = np.median(orgec_to_phopt[orgec])
        else:
            if len(org_to_phopt[org])>0:
                ph = np.median(org_to_phopt[org])
            else:
                ph = None
    tempcol.append(temp)
    phcol.append(ph)

df['temperature'] = tempcol
df['ph'] = phcol

def _split_reaction(reaction):
    reacs, prods = reaction.split(" >> ")
    reacs = reacs.split(' + ')
    prods = prods.split(' + ')
    return reacs, prods

rxnsmi_col = []
mwdiff_col = []
mw_col = []
all_rxnsmis = []
for substrate, subsmi, sub_to_reactions in tqdm(zip(df.substrate,
                                                    df.substrate_smiles,
                                                    df.reactions)):
    rxnsmis = []
    for sub, reactions in sub_to_reactions.items():
        try:
            subsminow = metabolite_inchi_smiles_dic.loc[sub].smiles
        except KeyError:
            continue
        if subsminow!=subsmi: 
            continue 
        for reaction in reactions:
            reacs, prods = _split_reaction(reaction)
            reacsmis = []
            prodsmis = []
            for reac in reacs:
                if reac in metabolite_inchi_smiles_dic.index: 
                    reacsmis.append(metabolite_inchi_smiles_dic.loc[reac].smiles)
            for prod in prods:
                if prod in metabolite_inchi_smiles_dic.index: 
                    prodsmis.append(metabolite_inchi_smiles_dic.loc[prod].smiles)
                    
            try:
                reac_mw = sum_mwts(reacsmis)
            except:
                reac_mw = None
                
            try:
                prod_mw = sum_mwts(prodsmis)
            except:
                prod_mw = None

            reacsmis_ = []
            for smi in reacsmis:
                if pd.isna(smi): 
                    continue
                else:
                    reacsmis_.append(smi)
            prodsmis_ = []
            for smi in prodsmis:
                if pd.isna(smi): 
                    continue
                else:
                    prodsmis_.append(smi)
                    
            reaction_smiles = f'{".".join(reacsmis_)}>>{".".join(prodsmis_)}'
            
            if pd.isna(reac_mw) or pd.isna(prod_mw): 
                rxnsmis.append((reaction_smiles,100000,0))
            else:
                rxnsmis.append((reaction_smiles,abs(reac_mw-prod_mw),reac_mw+prod_mw))
    
    rxn_now = None
    mwdiff_now = None
    mw_now = None
    # if something found, add 
    if len(rxnsmis)>0:
        rxnsmis_sorted = sorted(rxnsmis,key=sort_by_second)
        for rxnsmi, mwdiff,mw in rxnsmis_sorted:
            reacside, prodside = rxnsmi.split('>>')
            if subsmi in reacside: 
                rxn_now = rxnsmi
                mwdiff_now = mwdiff
                mw_now = mw
                break
                               
    mwdiff_col.append(mwdiff_now)
    mw_col.append(mw_now)
    rxnsmi_col.append(rxn_now)
    all_rxnsmis.append(rxnsmis)
    
df['reaction_smiles'] = rxnsmi_col
df['reaction_mw_difference'] = mwdiff_col
df['reaction_mw'] = mw_col

unicol = []
for ind, row in tqdm(df.iterrows()):
    if not type(row.uniprot) is list:
        if pd.isna(row.uniprot):
            unicol.append(None)
            continue
    unis = np.array(row.uniprot).flatten()
    try:
        unicol.append(';'.join(unis))
    except TypeError:
        unis_ = []
        for each in unis:
            unis_.extend(each)
        unis_ = np.array(unis_).flatten()
        unicol.append(';'.join(unis_))
            
df['uniprot'] = unicol

seqcol = []
uniset = set()
for uni in df.uniprot:
    if pd.isna(uni): continue
    else: 
        unis = uni.split(';')
        for each in unis:
            uniset.add(each)
            
def _get_sequence(uni): 
    r = requests.get(f'https://rest.uniprot.org/uniprotkb/{uni}.fasta')
    if r.status_code==200:
        lines = r.text.split('\n')
        seq = ''.join(lines[1:])
        return seq
    
outputs = Parallel(n_jobs=30, verbose=5)(delayed(_get_sequence)(uni) for uni in uniset)
uni_to_seq = {}
for uni, seq in zip(uniset, outputs):
    uni_to_seq[uni] = seq
    
seqrow = []
seq_srcrow = []
for ind, row in tqdm(df.iterrows()):
    uni = row.uniprot
    if pd.isna(uni): 
        seqrow.append(None)
        seq_srcrow.append(None)
    else: 
        unis = row.uniprot.split(';')
        seqs = []
        for uni in unis:
            if uni in uni_to_seq: 
                seqs.append(uni_to_seq[uni])
            else: 
                continue
        set_seqs = set(seqs)
        if len(set_seqs)>=1:
            seqrow.append(';'.join(set_seqs))
            seq_srcrow.append('brenda')
        else:
            seqrow.append(None)
            seq_srcrow.append(None)
            
df['sequence'] = seqrow
df['sequence_source'] = seq_srcrow
noseq_pairs = []
for ind, row in df.iterrows():
    if pd.isna(seq) or seq.strip()=='':
        if not pd.isna(tax) and not pd.isna(ec):
            seq = row.sequence
            tax = str(int(row.taxonomy_id))
            ec = row.ec
            noseq_pairs.append(tax+'__'+ec)
            
noseq_pairs = set(noseq_pairs)
noseq_tax_ec = []
for each in noseq_pairs:
    tax, ec = each.split('__')
    noseq_tax_ec.append((tax, ec))
    
def get_url(taxid,ec):
    return f"https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Csequence&format=tsv&query=%28%28organism_id%3A{taxid}%29+AND+%28ec%3A{ec}%29%29"

import requests
def fetch_and_process(items):
    ec,tax = items
    url = get_url(ec, tax)
    key = f'{tax}__{ec}.tsv'
    keydir = f'{args.processed_dir}/tax_ec_seqdata/'
    failed = False
    dic = {}
    if os.path.exists(keydir+key): 
         _df = pd.read_csv(keydir+key, sep='\t')
    else:
        r = requests.get(url)
        if r.status_code==200:
            text = r.text
            try:
                f = open(keydir + key, 'w')
                f.write(text)
                f.close()
            except:
                pass
            _df = pd.read_csv(keydir+key, sep='\t')
        else:
            failed = True
    if not failed:
        for entry, seq in zip(_df.Entry, _df.Sequence):
            dic[entry] = seq
        
    return dic
    
noseq_tax_ec = list(noseq_tax_ec)
seq_dicts = Parallel(n_jobs=30, verbose=100)(delayed(fetch_and_process)(items) for items in noseq_tax_ec)

tax_ec_seqs_dict = {}

for (tax,ec), seqs in zip(noseq_tax_ec, seq_dicts):
    tax_ec_seqs_dict[str(tax)+'__'+ec] = seqs
    
f = open(f'{args.processed_dir}/tax_ec_seqs_dict.json','w')
f.write(json.dumps(tax_ec_seqs_dict,indent=True))
f.close()

seqrow = []
srcrow = []
unirow = []
for ind, row in tqdm(df.iterrows()):
    seq = row.sequence
    src = row.sequence_source
    uni = row.uniprot
    if pd.isna(seq) and not pd.isna(row.taxonomy_id):
        tax = str(int(row.taxonomy_id))
        ec = str(row.ec)
        pair = tax+'__'+ec
        if pair in tax_ec_seqs_dict:
            seq_dict = tax_ec_seqs_dict[pair]
        else:
            seq_dict = {}
        if len(seq_dict)>0: 
            unis = []
            seqs = []
            for uni, seq in seq_dict.items():
                unis.append(uni)
                seqs.append(seq)
            seqrow.append(';'.join(seqs))
            unirow.append(';'.join(unis))
            srcrow.append('uniprot_search')
        else: 
            seqrow.append(seq)
            unirow.append(uni)
            srcrow.append('brenda')
    else:
        seqrow.append(seq)
        unirow.append(uni)
        srcrow.append('brenda')

df['uniprot'] = unirow
df['sequence'] = seqrow
df['sequence_source'] = srcrow

nseq_col = []
for ind, row in tqdm(df.iterrows()):
    if pd.isna(row.sequence): 
        nseq_col.append(0)
        seqcol.append(None)
    else:
        seqs = row.sequence.split(';')
        nseq_col.append(len(seqs))
    
df['n_sequence'] = nseq_col
        
df.dropna(subset=['sequence'],inplace=True)
df.reset_index(inplace=True, drop=True)

dfseq1 = df[df.n_sequence==1]
dfseq1.drop(columns=['reactions'],inplace=True)
dfseq1.reset_index(inplace=True, drop=True)
dfseq1.to_csv(f'{args.processed_dir}/brenda_processed_all_singleSeqs.csv')

dfseq1_wt = dfseq1[dfseq1.enzyme_type=='wild type']
dfseq1_wt.reset_index(inplace=True, drop=True)
dfseq1_wt.to_csv(f'{args.processed_dir}/brenda_processed_wt_singleSeqs.csv')

dfseq2 = df[df.n_sequence>=1]
dfseq2 = df[df.n_sequence<=10]
dfseq2.drop(columns=['reactions'],inplace=True)
dfseq2.reset_index(inplace=True, drop=True)
dfseq2.to_csv(f'{args.processed_dir}/brenda_processed_all_multipleSeqs.csv')

dfseq2_wt = dfseq2[dfseq2.enzyme_type=='wild type']
dfseq2_wt.reset_index(inplace=True, drop=True)
dfseq2_wt.to_csv(f'{args.processed_dir}/brenda_processed_wt_multipleSeqs.csv')
