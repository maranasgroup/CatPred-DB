This repository contains the datasets and scripts pertaining to the publication "CatPred: A comprehensive framework for deep learning in vitro enzyme kinetic parameters kcat, Km and Ki"
[![DOI](https://img.shields.io/badge/DOI-10.1101/2024.03.10.584340-blue)](https://www.biorxiv.org/content/10.1101/2024.03.10.584340v2)

<details><summary><b>Citation</b></summary>
CatPred biorxiv pre-print:
	
```bibtex
@article {Boorla2024.03.10.584340,
	author = {Veda Sheersh Boorla and Costas D. Maranas},
	title = {CatPred: A comprehensive framework for deep learning in vitro enzyme kinetic parameters kcat, Km and Ki},
	elocation-id = {2024.03.10.584340},
	year = {2024},
	doi = {10.1101/2024.03.10.584340},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/03/26/2024.03.10.584340},
	eprint = {https://www.biorxiv.org/content/early/2024/03/26/2024.03.10.584340.full.pdf},
	journal = {bioRxiv}
}
```
</details>

<details open><summary><b>Table of contents</b></summary>
	
- [CatPred-DB datasets](#datasets)
- [CatPred-DB processing pipeline](#pipeline)
- [License](#license)
</details>

## CatPred-DB datasets <a name="datasets"></a>

### Directory information for datasets

├── datasets
    ├── processed               		# processed datasets one each for kcat, Km and Ki
    ├── splits                  		# training/validation/test dataset splits for kcat, Km and Ki
    ├── all_natural_metabolite_names.json 		# compiled names of natural metabolites from BRENDA
    ├── metabolite_inchi_smiles_brenda_pubchem.tsv 	# compiled list of inchi and smiles strings for brenda molecules
    

# Obtaining raw data
1. BRENDA raw data has been obtained from their website. Recently, they made .json format download for the database available. In this work, 'brenda_2022_2.json' is used. Also, raw data for compounds was obtained from https://www.brenda-enzymes.org/search_result.php?a=13 and placing a blank query. Thanks to samgoldman97 for the trick. Even with this some brenda compounds didn't have a known SMILES/InChi string. So, we also used the Pubchem id exchange service to map names (synonyms) to SMILES. All these are saved in ./data/raw/brenda/
2. SABIO-rk raw data has been obtained using sbml exports from their website. Because only a 100 entries are shown in their website at a time, we manually exported a lot of sbml files from their website. Also, we scraped individual entry html files. These form the raw data. Saved in ./data/raw/sabio/

# Processing raw data
1. Processing BRENDA is relatively straightforward because everything is in the json file except SMILES strings of substrates and Sequences of proteins. See scripts/data/processing/ for details.
2. SABIO sbml files are a bit more tricky to process. We used an sbml parser from libsbml. Also, BeautifulSoup to parse html files of individual entries. See scripts/data/processing/ for details.

# Merging processed data
1. Since data comes from two different sources, they need to be merged to make sure there are no unexpected duplicates. Also, each 'sequence,smiles' pair can appear in several entries (possibly with different parameter values). So, these have to be handled appropriately. See scripts/data/merging/ for details.

# Splitting merged data
1. The final merged data is split for training/validation/testing appropriately for each parameter separately. See scripts/data/splitting/
