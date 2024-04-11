import argparse
import re
from bs4 import BeautifulSoup
import sys
import os.path
from libsbml import *
import json
from tqdm import tqdm
from rdkit.Chem import Descriptors
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import json
import sys
import os
import os.path

import warnings
warnings.filterwarnings("ignore")

from libsbml import *
import json
import numpy as np
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from rdkit import Chem
from joblib import Parallel, delayed
import ipdb
import ete3

def get_inchi_smiles(html_string):

    soup = BeautifulSoup(html_string, 'lxml') # Parse the HTML as a string

    table = soup.find_all('table')[0] # Grab the first table

    tds = table.find_all('td')

    items_all = []
    for td in tds:
        items = td.find_all('span')
        if not items: continue
        else: items_all.extend(items)

    inchis = []
    smiles = []
    for item in items_all:
        val = item.attrs['id']
        if 'Inchi' in val:
            inchis.append(item.text)
        elif 'Smiles' in val: 
            smiles.append(item.text)

    return inchis, smiles

def has_regulator(html_string):
    soup = BeautifulSoup(html_string, 'lxml') # Parse the HTML as a string

    name_to_cid = {}
    
    try:
        table = soup.find_all('table')[0] # Grab the first table
    except:
        return name_to_cid
    
    tds = table.find_all('td')
    
    for td in tds:
        if 'Modifier-Inhibitor' in td or 'Modifier-Activator' in td:
            return True
    
    return False
    
def get_name_to_cid(html_string):
    soup = BeautifulSoup(html_string, 'lxml') # Parse the HTML as a string

    name_to_cid = {}
    
    try:
        table = soup.find_all('table')[0] # Grab the first table
    except:
        return name_to_cid
    
    tds = table.find_all('td')
    
    for td in tds:
        items = td.find_all('a')
        for item in items:
            if not 'compdetails.jsp?' in str(item): continue
            cid = item.attrs['href'].split('./compdetails.jsp?cid=')[-1]
            name = item.text.strip()
            name_to_cid[name] = cid

    return name_to_cid

def get_protein_units(html_string):
    soup = BeautifulSoup(html_string, 'lxml') # Parse the HTML as a string

    table = soup.find_all('table')[0] # Grab the first table

    tds = table.find_all('td')

    ourlines = []
    for td in tds:
        std = str(td)
        if not 'proteindetails.jsp?' in std: continue
        lines = std.split()
        for line in lines:
            if 'http://sabiork.h-its.org/proteindetails.jsp?' in line:
                ourlines.append(line)

    ourlines = list(set(ourlines))

    uniprots = []
    repeats = []

    for line in ourlines:
        try:
            if ')*' in line: item = re.findall(r'>[A-Z,0-9]+</a>[\)]*[\*,0-9]*',line)[0]
            else: item = re.findall(r'>[A-Z,0-9]+</a>[\)]*',line)[0]
        except:
            print(line)
        uni, extra = item.split('</a>')
        uni = uni.split('>')[-1]
        extra = extra.split(';')[0]
        if '*' in extra:
            repeat = extra.split('*')[-1]
        else:
            repeat = '1'
        uniprots.append(uni)
        repeats.append(repeat)

    return uniprots, repeats


def parse(filename):
    document = readSBML(filename)

    errors = document.getNumErrors()

    m = document.getModel()
 
    species_data = {
    }
    # map species id to name
    # Species
    #print("----- " + 'Species' + "----- ")
    for i in range(0, m.getNumSpecies()):
        data = {
        }
        sp = m.getSpecies(i)
        sp_str = str(sp)[1:-1] # 1, -1 to remove < >
        items = sp_str.split('Species')[1].strip().split()
        sp_id = items[0]
        sp_name = ' '.join(items[1:])
        #print(sp_id, ':', sp_name[1:-1])
        data['name'] = sp_name[1:-1]
        annot_lines = sp.annotation_string.split('\n')
        # #print(annot_lines)
        chebi_ids = []
        kegg_ids = []
        uni_ids = []
        for line in annot_lines:
            if 'rdf:resource' in line:
                link = line.split('rdf:resource')[1][2:-3]
                some_id = link.split('/')[-1]
                if 'chebi' in link:
                    chebi_ids.append(some_id)
                elif 'kegg' in link:
                    kegg_ids.append(some_id)
                elif 'uniprot' in link:
                    uni_ids.append(some_id)
                #print(link)
        data['chebi'] = chebi_ids
        data['kegg'] = kegg_ids
        data['uniprot'] = uni_ids
    
        species_data[sp_id] = data
 
   # FunctionDefinition
    # for i in range (0, m.getNumFunctionDefinitions()):
    #     sp = m.getFunctionDefinition(i)
    #     #printAnnotation(sp)
 
    unit_data = {} # map unit id to name
    # UnitDefinition
    #print("----- " + "unitDefinition" + "----- ")
    for i in range (0, m.getNumUnitDefinitions()):
        sp = m.getUnitDefinition(i)
        unit_str = str(sp)[1:-1] # 1, -1 to remove < >
        unitid, unitname = unit_str.split('UnitDefinition')[1].strip().split()
        unit_data[unitid] = unitname[1:-1] # 1, -1 to remove " "
        #print(unitid, ':', unitname[1:-1])

    rxn_data = []
    for i in range(0, m.getNumReactions()):
        data = {}
        re = m.getReaction(i)
        annot_lines = re.annotation_string.split('\n')
        re_links = []
        #print("----- " + 'Reaction' + ' id: ' + re.id + " ----- ")
        ec_ids = []
        tax_ids = []
        sabio_rxn_ids = []
        for line in annot_lines:
            if 'rdf:resource' in line:
                link = line.split('rdf:resource')[1][2:-3]
                re_links.append(link)
                some_id = link.split('/')[-1]
                if 'ec' in link:
                    ec_ids.append(some_id)
                elif 'taxonomy' in link:
                    tax_ids.append(some_id)
                elif 'sabiork' in link:
                    sabio_rxn_ids.append(some_id)
                #print(link)
        data['ec'] = ec_ids
        data['taxonomy'] = tax_ids
        data['sabio_reaction'] = sabio_rxn_ids

        # SpeciesReference (Reactant)

        #print("----- " + 'Reactants' + " ----- ")
        reactant_species_ids = []
        reactant_stoichs = []
        for j in range(0, re.getNumReactants()):
            rt = re.getReactant(j)
            reactant_stoichs.append(rt.getStoichiometry())
            reactant_species_ids.append(rt.getSpecies())
            #print(rt.getStoichiometry(), rt.getSpecies())
 
        data['reactants'] = {'species': reactant_species_ids,'stoichiometry':reactant_stoichs}
        # SpeciesReference (Product)
 
        #print("----- " + 'Products' + " ----- ")
        product_species_ids = []
        product_stoichs = []
        for j in range(0, re.getNumProducts()):
            rt = re.getProduct(j)
            product_stoichs.append(rt.getStoichiometry())
            product_species_ids.append(rt.getSpecies())
            #print(rt.getStoichiometry(), rt.getSpecies())
 
        data['products'] = {'species': product_species_ids,'stoichiometry':product_stoichs}

        # ModifierSpeciesReference (Modifiers)
 
        #print("----- " + 'Modifiers' + " ----- ")
        modifier_species_ids = []
        for j in range(0, re.getNumModifiers()):
            md = re.getModifier(j)
            modifier_species_ids.append(md.getSpecies())
            #print(md.getSpecies())

        data['modifiers'] = modifier_species_ids

        # KineticLaw
 
        if re.isSetKineticLaw():
            kl = re.getKineticLaw()
            #print("----- " + 'Kinetic Law' + " ----- ")
            annot_lines = kl.annotation_string.split('\n')
            kl_data = {}
            temperature = ''
            temperature_unit = ''
            ph = ''
            buffer = ''
            for line in annot_lines:
                tag = 'sbrk:startValueTemperature'
                if tag in line:
                    items = line.split(tag)
                    temperature = items[1][1:-2]
                tag = 'sbrk:temperatureUnit'
                if tag in line:
                    items = line.split(tag)
                    temperature_unit = items[1][2:-2]
                tag = 'sbrk:startValuepH'
                if tag in line:
                    items = line.split(tag)
                    ph = items[1][1:-2]
                tag = 'sbrk:buffer'
                if tag in line:
                    items = line.split(tag)
                    buffer = items[1][1:-2]
                if 'rdf:resource' in line:
                    link = line.split('rdf:resource')[1][2:-3]
                    some_id = link.split('/')[-1]
                    if 'pubmed' in link:
                        kl_data['pubmed'] = some_id
                    elif 'sabiork' in link:
                        kl_data['sabio_kinetic_law'] = some_id

            #print("----- " + 'Experimental conditions' + " ----- ")
            #print('Temperature:', temperature, temperature_unit)
            #print('pH:', ph)
            #print('Buffer:', buffer)
            kl_data['temperature'] = temperature
            kl_data['ph'] = ph
            kl_data['buffer'] = buffer

            # Parameter
            parameter_ids = []
            parameter_values = []
            parameter_units = []

            #print("----- " + 'Parameters' + " ----- ")
            for j in range(0, kl.getNumParameters()):
                pa = kl.getParameter(j)
                parameter_ids.append(pa.id)
                parameter_values.append(pa.getValue())
                parameter_units.append(pa.getUnits())
                #print(pa.id, ':', pa.getValue(), pa.getUnits())

            kl_data['parameters'] = {
                'id': parameter_ids,
                'values': parameter_values,
                'units': parameter_units
            }

        data['kinetic_law'] = kl_data
        rxn_data.append(data)

    return species_data, unit_data, rxn_data 

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rawdata_dir", help="directory where raw data is stored", type=str, required=True)
    parser.add_argument("--processed_dir", help="directory where to store processed data", type=str, required=True)
    parser.add_argument("--out_file_name", help="out file name", 
                        type=str, required=True)

    args, unparsed = parser.parse_known_args()
    parser = argparse.ArgumentParser()

    return args

args = parse_args()
N_PROCS = 30

sbml_data_dir = os.path.abspath(args.rawdata_dir)
processed_dir = os.path.abspath(args.processed_dir)
out_path = processed_dir + args.out_file_name

# these xml files were downloaded using SABIO-RK web search interface for wild-type enzymes
         
sbml_files = [os.path.join(sbml_data_dir, file) for file in os.listdir(sbml_data_dir) if file.endswith('.xml')]
species_data_all = {}
unit_data_all = {}
rxn_data_all = []
for sbml_file in tqdm(sbml_files):
    species_data, unit_data, rxn_data = parse(sbml_file.strip())
    species_data_all.update(species_data)
    unit_data_all.update(unit_data)
    rxn_data_all.extend(rxn_data)

f = open(os.path.join(processed_dir, 'sabio_sbml_species_dict.json'),'w')
f.write(json.dumps(species_data_all,indent=True))
f.close()
f = open(os.path.join(processed_dir, 'sabio_sbml_unit_dict.json'),'w')
f.write(json.dumps(unit_data_all,indent=True))
f.close()
f = open(os.path.join(processed_dir, 'sabio_sbml_reaction_dict.json'),'w')
f.write(json.dumps(rxn_data_all,indent=True))
f.close()

rxn_data = json.load(open(os.path.join(processed_dir, 'sabio_sbml_reaction_dict.json')))
species_data = json.load(open(os.path.join(processed_dir, 'sabio_sbml_species_dict.json')))
unit_data = json.load(open(os.path.join(processed_dir, 'sabio_sbml_unit_dict.json')))

kl_sabio_id_col = []
pubmed_col = []
uni_col = []
reactant_col = []
reactant_stoich_col = []
product_col = []
product_stoich_col = []

param_species_col = []

param_col = []
param_value_col = []
param_unit_col = []
temperature_col = []
ph_col = []
buffer_col = []
ec_col = []
tax_col = []
rxn_sabio_id_col = []

for rxn in rxn_data:
    for param,value,unit in zip(rxn['kinetic_law']['parameters']['id'],rxn['kinetic_law']['parameters']['values'],rxn['kinetic_law']['parameters']['units']):
        param_col.append(param)
        param_value_col.append(value)
        param_unit_col.append(unit)
        rxn_sabio_id_col.append(';'.join(rxn['sabio_reaction']))
        tax_col.append(';'.join(rxn['taxonomy']))
        ec_col.append(';'.join(rxn['ec']))
        uni = False
        unis = []
        for mod in rxn['modifiers']:
            if mod.startswith('ENZ') and mod in species_data:
                if len(species_data[mod]['uniprot'])>0:
                    uni = True
                    unis.append(';'.join(species_data[mod]['uniprot']))
                
        if len(unis)==0: 
            uni_col.append(None)
        else:
            uni_col.append(';'.join(unis))

        if 'pubmed' in rxn['kinetic_law']:
            pubmed_col.append(rxn['kinetic_law']['pubmed'])
        else:
            pubmed_col.append(None)

        if 'sabio_kinetic_law' in rxn['kinetic_law']:
            kl_sabio_id_col.append(rxn['kinetic_law']['sabio_kinetic_law'])
        else:
            kl_sabio_id_col.append(None)

        reactant_col.append(rxn['reactants']['species'])
        reactant_stoich_col.append(rxn['reactants']['stoichiometry'])

        product_col.append(rxn['products']['species'])
        product_stoich_col.append(rxn['products']['stoichiometry'])

        buffer_col.append(rxn['kinetic_law']['buffer'])
        temperature_col.append(rxn['kinetic_law']['temperature'])
        ph_col.append(rxn['kinetic_law']['ph'])

sabio_df = pd.DataFrame()
sabio_df['uniprot'] = uni_col
sabio_df['ec'] = ec_col
sabio_df['taxonomy_id'] = tax_col
sabio_df['param'] = param_col
sabio_df['param_value'] = param_value_col
sabio_df['param_unit'] = param_unit_col
sabio_df['reactant_species'] = reactant_col
sabio_df['product_species'] = product_col
sabio_df['reactant_stoichiometry'] = reactant_stoich_col
sabio_df['product_stoichiometry'] = product_stoich_col

sabio_df['temperature'] = temperature_col
sabio_df['ph'] = ph_col
sabio_df['buffer'] = buffer_col
sabio_df['pubmed_id'] = pubmed_col
sabio_df['sabio_reaction_id'] = rxn_sabio_id_col
sabio_df['sabio_kinetic_law_id'] = kl_sabio_id_col

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

ec_words = []
for ind, row in sabio_df.iterrows():
    words = get_ec_words(row.ec)
    ec_words.append(words)
for col in ec_embed_cols:
    col_values = [ec_words[i][col] for i in range(len(sabio_df))]
    sabio_df[col] = col_values

tax_words = []
for ind, row in sabio_df.iterrows():
    words = get_tax_words(row.taxonomy_id, ncbi)
    tax_words.append(words)
for col in tax_embed_cols:
    col_values = []
    for i in range(len(sabio_df)):
        if col in tax_words[i]:
            col_values.append(tax_words[i][col])
        else:
            col_values.append('UNK')
    sabio_df[col] = col_values
    
# ipdb.set_trace()
# Add species names

def convert_unit_kcat(kcat,unit):
    if unit=='minwedgeone': 
        kcat = kcat*60
    elif unit=='swedgeone': 
        kcat = kcat
    elif unit=='hwedgeone': 
        kcat = kcat*60*60
    else:
        kcat = None
    return kcat

def convert_unit_km(km,unit):
    if unit=='M': 
        km = km*1000
    elif unit=='microM': 
        km = km*0.001
    elif unit=='nM': 
        km = km*0.001*0.001
    elif unit=='mM':
        km = km
    else:
        km = km
    return km

# def convert_unit_kcat_km(kcat_km, unit):
#     if unit=='mMwedgeoneswedgeone':
#         kcat_km = kcat_km
#     elif unit=='Mwedgeoneminwedgeone':
#         kcat_km = kcat_km/60.
#         kcat_km = kcat_km/
#     ['','Mwedgeoneswedgeone','mMwedgeoneminwedgeone','',
#                           'microMwedgeoneminwedgeone','nMwedgeoneminwedgeone','nMwedgeoneswedgeone']

reactant_names_col = []
product_names_col = []
param_names_col = []
param_values_col = []
units = []
param_species_col = []
param_species_name_col = []
total = 0
u = 0
for reactants,products,param,param_unit,param_value in zip(
            sabio_df.reactant_species,
            sabio_df.product_species,
            sabio_df.param,
            sabio_df.param_unit,
            sabio_df.param_value):
    
    if param=='kcat':
        par = 'kcat'
        val = convert_unit_kcat(param_value, param_unit)
    elif param.startswith('Km_'):
        _, species = param.split('Km_')
        par = 'Km'
        val = convert_unit_km(param_value, param_unit)
    elif param.startswith('Ki_'):
        _, species = param.split('Ki_')
        par = 'Ki'
        val = convert_unit_km(param_value, param_unit)
    elif param.startswith('IC50_'):
        _, species = param.split('IC50_')
        par = 'IC50'
        val = convert_unit_km(param_value, param_unit)
    elif param.startswith('kcat_Km_'):
        _, species = param.split('kcat_Km_')
        par = 'kcat_Km'
        units.append(param_unit)
        # total+=1
        # if param_unit in ['Mwedgeoneminwedgeone','Mwedgeoneswedgeone','mMwedgeoneminwedgeone','mMwedgeoneswedgeone',
        #                   'microMwedgeoneminwedgeone','nMwedgeoneminwedgeone','nMwedgeoneswedgeone']: u+=1
        val = param_value#convert_unit_km(param_value, param_unit)
    else:
        par = None
        val = None
    param_values_col.append(val)
    param_species_col.append(species)
    param_names_col.append(par)
    reactant_names = []
    product_names = []
    for reactant in reactants:
        if pd.isna(reactant): 
            reactant_names.append(None)
            continue
        elif reactant in species_data:
            name = species_data[reactant]['name']
            reactant_names.append(name)
    for product in products:
        if pd.isna(product): 
            product_names.append(None)
            continue
        elif product in species_data:
            name = species_data[product]['name']
            product_names.append(name)
    if species in species_data:
        name = species_data[species]['name']
        param_species_name_col.append(name)
    else:
        param_species_name_col.append(None)
    reactant_names_col.append(reactant_names)
    product_names_col.append(product_names)

sabio_df['reactant_names'] = reactant_names_col
sabio_df['product_names'] = product_names_col
sabio_df['param_name'] = param_names_col
sabio_df['param_value_stdunit'] = param_values_col
sabio_df['param_species'] = param_species_col
sabio_df['param_species_name'] = param_species_name_col

print(sabio_df.columns)

from joblib import Parallel, delayed
import requests

def download_sabio_kl(kl):
    r = requests.get(f'https://sabiork.h-its.org/kindatadirectiframe.jsp?kinlawid={kl}')
    if r.status_code==200: 
        html_string = r.text
        f = open(f'{sbml_data_dir}/kl_html/{kl}.html','w')
        f.write(html_string)
        f.close()
        return True
    else:
        return False
    
kls_to_download = set()
for klid in tqdm(sabio_df.sabio_kinetic_law_id.unique()):
    html_file = f'{sbml_data_dir}/kl_html/{klid}.html'
    if not os.path.exists(html_file): 
        kls_to_download.add(klid)
        
outputs = Parallel(n_jobs=30, verbose=5)(delayed(download_sabio_kl)(kl) for kl in kls_to_download)

def download_sabio_cid(cid):
    r = requests.get(f'https://sabiork.h-its.org/compdetails.jsp?cid={cid}')
    if r.status_code==200: 
        html_string = r.text
        f = open(f'{sbml_data_dir}/cid_html/{cid}.html','w')
        f.write(html_string)
        f.close()
        return True
    else:
        return False
    os.path.exists(html_file)
    
klids_regulated = []
name_to_cid = {}
for klid in tqdm(sabio_df.sabio_kinetic_law_id.unique()):
    html_file = f'{sbml_data_dir}/kl_html/{klid}.html'
    try:
        with open(html_file) as f:
            html_text = f.read()
            name_to_cid.update(get_name_to_cid(html_text))
            if has_regulator(html_text): klids_regulated.append(klid)
    except:
        continue

cids_to_download = set()
for name, cid in name_to_cid.items():
    html_file = f'{sbml_data_dir}/cid_html/{cid}.html'
    if not os.path.exists(html_file): 
        cids_to_download.add(cid)

outputs = Parallel(n_jobs=30, verbose=5)(delayed(download_sabio_cid)(cid) for cid in cids_to_download)

for klid in tqdm(sabio_df.sabio_kinetic_law_id.unique()):
    html_file = f'{sbml_data_dir}/kl_html/{klid}.html'
    try:
        with open(html_file) as f:
            html_text = f.read()
            name_to_cid.update(get_name_to_cid(html_text))
            if has_regulator(html_text): klids_regulated.append(klid)
    except:
        continue

regulated_col = []
for klid in sabio_df.sabio_kinetic_law_id:
    if klid in klids_regulated: regulated_col.append(True)
    else: regulated_col.append(None)
    
sabio_df['has_regulator'] = regulated_col

all_metabolite_names = set()
name_to_smiles = {}
name_to_inchi = {}
nosmi_names = set()

for reactants, products,species in zip(sabio_df.reactant_names,sabio_df.product_names,sabio_df.param_species_name):
    for each in reactants+products:
        if not pd.isna(each): all_metabolite_names.add(each)
    if not pd.isna(species): all_metabolite_names.add(species)

def inchi_to_smiles(inchi): 
    smi = Chem.MolToSmiles(Chem.MolFromInchi(inchi))
    return smi if not smi is None else 'None'

for name in tqdm(all_metabolite_names):
    if name in name_to_cid:
        cid = name_to_cid[name]
    else:
        nosmi_names.add(name)
        continue
    html_file = f'{sbml_data_dir}/cid_html/{cid}.html'
    with open(html_file) as f:
        html_text = f.read()
        inchi, smiles = get_inchi_smiles(html_text)
        if smiles: name_to_smiles[name] = smiles[0]
        if inchi: name_to_inchi[name] = inchi[0]

for name, inchi in name_to_inchi.items():
    if name in name_to_smiles: continue
    name_to_smiles[name] = inchi_to_smiles(inchi)

name_to_smiles_can = {}
for name,smi in tqdm(name_to_smiles.items()):
    if pd.isna(smi):
        continue
    else:
        try:
            name_to_smiles_can[name] = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        except:
            name_to_smiles_can[name] = None

name_to_smiles = name_to_smiles_can.copy()
metabolite_inchi_smiles_dic = pd.read_csv(f'{sbml_data_dir}/metabolite_inchi_smiles_brenda_pubchem.tsv',sep='\t',
                                          index_col='metabolite')
for name,smi in tqdm(name_to_smiles.items()):
    if pd.isna(smi):
        if name in metabolite_inchi_smiles_dic.index:
            name_to_smiles[name] = metabolite_inchi_smiles_dic.loc[name].smiles

reac_smiles_col = []
prod_smiles_col = []
param_species_smiles_col = []

for reactants, products,param_species in zip(sabio_df.reactant_names,sabio_df.product_names,sabio_df.param_species_name):
    reac_smiles = []
    prod_smiles = []
    for reac in reactants:
        if pd.isna(reac):
            smi = None
        elif not reac in name_to_smiles:
            nosmi_names.add(reac)
            smi = None
        else:
            smi = name_to_smiles[reac]
        reac_smiles.append(smi)
    reac_smiles_col.append(reac_smiles)
    for prod in products:
        if pd.isna(prod):
            smi = None
        elif not prod in name_to_smiles:
            nosmi_names.add(prod)
            smi = None
        else:
            smi = name_to_smiles[prod]
        prod_smiles.append(smi)
    prod_smiles_col.append(prod_smiles)
    if param_species in name_to_smiles:
        smi = name_to_smiles[param_species]
    else:
        nosmi_names.add(param_species)
        smi = None
    param_species_smiles_col.append(smi)
    
sabio_df['reactant_smiles'] = reac_smiles_col
sabio_df['product_smiles'] = prod_smiles_col
sabio_df['param_species_smiles'] = param_species_smiles_col
print(sabio_df.columns)

f = open(f'{sbml_data_dir}/pubchem_nosmiles.txt','w')
for sub in nosmi_names:
    if not pd.isna(sub):f.write(sub+'\n')
f.close()

# get substrate to smiles mapping from https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi

sub_to_smiles_data = open(f'{sbml_data_dir}/pubchem_sub_to_smiles.txt').read()
lines = sub_to_smiles_data.split('\n')
for line in lines:
    if not line.strip(): continue
    sub, smi = line.split('\t')
    if smi: name_to_smiles[sub] = smi

smi_to_mw_dic = {}
for smi in tqdm(name_to_smiles.values()):
    mw = None
    try:
        mol = Chem.MolFromSmiles(smi)
        mw = Descriptors.MolWt(mol)
    except:
        pass
    smi_to_mw_dic[smi] = mw
    
def sum_mwts(smis):
    mw = 0
    for smi in smis:
        mw+=smi_to_mw_dic[smi]
    return mw

mwdiff_col = []
mw_col = []
reaction_smiles_col = []
param_species_smiles_col = []

for reactants, products,param_species, rst, pst in zip(sabio_df.reactant_names,sabio_df.product_names,sabio_df.param_species_name, sabio_df.reactant_stoichiometry, sabio_df.product_stoichiometry):
    rxnsmi_now = None
    mwdiff_now = None
    mw_now = None
    sminow = None
    reac_smiles = []
    prod_smiles = []
    
    if len(reactants)!=rst: 
        rst = [1]*len(reactants)
    if len(products)!=pst: 
        pst = [1]*len(products)
    
    for reac,r in zip(reactants,rst):
        if pd.isna(reac):
            continue
        elif not reac in name_to_smiles:
            continue
        elif not name_to_smiles[reac] is None:
            try:
                smi = '.'.join([name_to_smiles[reac]]*r)
            except:
                pass
        reac_smiles.append(smi)
        
    reac_smiles_col.append(reac_smiles)
    for prod,p in zip(products,pst):
        if pd.isna(prod):
            continue
        elif not prod in name_to_smiles:
            continue
        elif not name_to_smiles[prod] is None:
            try:
                smi = '.'.join([name_to_smiles[prod]]*p)
            except:
                continue
        prod_smiles.append(smi)
        
    prod_smiles_col.append(prod_smiles)
    try:
        reac_mw = sum_mwts(reac_smiles)
    except:
        reac_mw = None

    try:
        prod_mw = sum_mwts(prod_smiles)
    except:
        prod_mw = None
        
    reacsmis_ = []
    for smi in reac_smiles:
        if pd.isna(smi): 
            continue
        else:
            reacsmis_.append(smi)
    prodsmis_ = []
    for smi in prod_smiles:
        if pd.isna(smi): 
            continue
        else:
            prodsmis_.append(smi)

    rxnsmi_now = f'{".".join(reacsmis_)}>>{".".join(prodsmis_)}'
    
    if pd.isna(reac_mw) or pd.isna(prod_mw): 
        mwdiff_now = 100000
        mw_now = 0
    else:
        rxnsmi_now = f'{".".join(reac_smiles)}>>{".".join(prod_smiles)}'
        mwdiff_now = abs(reac_mw-prod_mw)
        mw_now = reac_mw + prod_mw
        
    if param_species in name_to_smiles:
        sminow = name_to_smiles[param_species]
        
    mwdiff_col.append(mwdiff_now)
    mw_col.append(mw_now)
    reaction_smiles_col.append(rxnsmi_now)
    param_species_smiles_col.append(sminow)
    
sabio_df['reaction_smiles'] = reaction_smiles_col
sabio_df['reaction_mw_difference'] = mwdiff_col
sabio_df['reaction_mw'] = mw_col
sabio_df['param_species_smiles'] = param_species_smiles_col
sabio_df = sabio_df[sabio_df.param_name.isin(['kcat','Km','Ki','kcat_Km'])]

sabio_df.drop(columns=['reactant_species','product_species','reactant_stoichiometry','product_stoichiometry','reactant_names','product_names','reactant_smiles','product_smiles'],inplace=True)

n_uniprot_col = []
for ind, row in sabio_df.iterrows():
    try:
        unis = row.uniprot.split(';')
        n_uniprot_col.append(len(unis))
    except:
        n_uniprot_col.append(0)

sabio_df['N_uniprot'] = n_uniprot_col

seqcol = []
uniset = set()
for uni in sabio_df.uniprot:
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
for ind, row in tqdm(sabio_df.iterrows()):
    uni = row.uniprot
    if pd.isna(uni): 
        seqrow.append(None)
    else: 
        unis = row.uniprot.split(';')
        seqs = []
        for uni in unis:
            if uni in uni_to_seq: seqs.append(uni_to_seq[uni])
            else: continue
        if len(seqs)>0:
            set_seqs = set()
            for s in seqs: 
                if not pd.isna(s): set_seqs.add(s)
            seqrow.append(';'.join(set_seqs))
        else:
            seqrow.append(None)
            
sabio_df['sequence'] = seqrow
sabio_df['sequence_source'] = 'sabio'

noseq_pairs = []
for ind, row in sabio_df.iterrows():
    seq = row.sequence
    if pd.isna(seq) or seq.strip()=='':
        noseq_pairs.append(row.taxonomy_id+'__'+row.ec)
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
    if not os.path.exists(keydir): os.mkdir(keydir)
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
for ind, row in tqdm(sabio_df.iterrows()):
    seq = row.sequence
    src = row.sequence_source
    uni = row.uniprot
    if pd.isna(seq) and not pd.isna(row.taxonomy_id):
        tax = str(row.taxonomy_id)
        ec = str(row.ec)
        pair = tax+'__'+ec
        if pair in tax_ec_seqs_dict:
            seq_dict = tax_ec_seqs_dict[pair]
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
            srcrow.append('sabio')
    else:
        seqrow.append(seq)
        unirow.append(uni)
        srcrow.append('sabio')

sabio_df['uniprot'] = unirow
sabio_df['sequence'] = seqrow
sabio_df['sequence_source'] = srcrow

nseq_col = []
for ind, row in tqdm(sabio_df.iterrows()):
    if pd.isna(row.sequence): 
        nseq_col.append(0)
        seqcol.append(None)
    else:
        seqs = row.sequence.split(';')
        nseq_col.append(len(seqs))
    
sabio_df['n_sequence'] = nseq_col

print(sabio_df.columns)

dfseq1 = sabio_df[sabio_df.n_sequence==1]
dfseq1.reset_index(inplace=True, drop=True)
dfseq1.to_csv(f'{args.processed_dir}/sabio_processed_wt_singleSeqs.csv')

dfseq2 = sabio_df[sabio_df.n_sequence>=1]
dfseq2.reset_index(inplace=True, drop=True)
dfseq2.to_csv(f'{args.processed_dir}/sabio_processed_wt_multipleSeqs.csv')

sabio_df.to_csv(out_path)
