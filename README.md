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
- [CatPred-DB scripts](#scripts)
- [License](#license)
</details>

## CatPred-DB datasets <a name="datasets"></a>

Datasets referenced in the publication can be found in these directories

    ./datasets/
    ├── processed               		# processed datasets one each for kcat, Km and Ki
    ├── splits                  		# training/validation/test dataset splits for kcat, Km and Ki
    ├── all_natural_metabolite_names.json 		# compiled names of natural metabolites from BRENDA
    ├── metabolite_inchi_smiles_brenda_pubchem.tsv 	# compiled list of inchi and smiles strings for brenda molecules
    
## CatPred-DB scripts <a name="datasets"></a>

Custom scripts used for processing, splitting and analysis can be found in these directories

    ./scripts/
    ├── processing               	# for processing brenda/sabio datasets and merging them
    ├── pdb                  		# for creating files required for datasets including structural information (pdb)
    ├── add_pdb_paths.py 		# adds pdb paths to a dataset file
    ├── count_natural.py 		# counts natural vs. non-natural entries
    ├── create_pdbrecords.py 		# create pdb records files to go along with datasets required for training
    ├── make_splits.py 			# create train/valid/test splits from processed datasets

## License <a name="license"></a>

This source code is licensed under the MIT license found in the `LICENSE` file
in the root directory of this source tree.
