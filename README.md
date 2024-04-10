Obtaining raw data

BRENDA raw data has been obtained from their website. Recently, they made .json format download for the database available. In this work, 'brenda_2022_2.json' is used. Also, raw data for compounds was obtained from https://www.brenda-enzymes.org/search_result.php?a=13 and placing a blank query. Thanks to samgoldman97 for the trick. Even with this some brenda compounds didn't have a known SMILES/InChi string. So, we also used the Pubchem id exchange service to map names (synonyms) to SMILES. All these are saved in ./data/raw/brenda/
SABIO-rk raw data has been obtained using sbml exports from their website. Because only a 100 entries are shown in their website at a time, we manually exported a lot of sbml files from their website. Also, we scraped individual entry html files. These form the raw data. Saved in ./data/raw/sabio/
Processing raw data

Processing BRENDA is relatively straightforward because everything is in the json file except SMILES strings of substrates and Sequences of proteins. See scripts/data/processing/ for details.
SABIO sbml files are a bit more tricky to process. We used an sbml parser from libsbml. Also, BeautifulSoup to parse html files of individual entries. See scripts/data/processing/ for details.
Merging processed data

Since data comes from two different sources, they need to be merged to make sure there are no unexpected duplicates. Also, each 'sequence,smiles' pair can appear in several entries (possibly with different parameter values). So, these have to be handled appropriately. See scripts/data/merging/ for details.
Splitting merged data

The final merged data is split for training/validation/testing appropriately for each parameter separately. See scripts/data/splitting/
