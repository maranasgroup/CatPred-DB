# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: MIT-0

"""
Parse PDB files to extract the coordinates of the 4 key atoms from AAs to
generate json records compatible to the LM-GVP model.

This script is intended for Fluorescence and Protease datasets from TAPE.
"""

import json
import os
import argparse
from collections import defaultdict
import pandas as pd
import numpy as np

from tqdm import tqdm
from joblib import Parallel, delayed
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import three_to_one

import xpdb
from contact_map_utils import parse_pdb_structure


def parse_args():
    """Prepare argument parser.

    Args:

    Return:

    """
    parser = argparse.ArgumentParser(
        description="Generate json records for GVP from pdb files"
    )
    parser.add_argument(
        "--data_file",
        help="Path to csv file containing pdbpaths",
        required=True,
    )
    parser.add_argument(
        "--target_col",
        help="target column in data",
        required=True,
    )
    parser.add_argument(
        "--smiles_col",
        help="smiles column in data",
        required=True,
    )
    parser.add_argument("--out_prefix", help="output prefix for json records")

    args = parser.parse_args()
    return args


def get_atom_coords(residue, target_atoms=["N", "CA", "C", "O"]):
    """Extract the coordinates of the target_atoms from an AA residue.

    Args:
        residue: a Bio.PDB.Residue object representing the residue.
        target_atoms: Target atoms which residues will be returned.

    Retruns:
        Array of residue's target atoms (in the same order as target atoms).
    """
    return np.asarray([residue[atom].coord for atom in target_atoms])


def structure_to_coords(struct, target_atoms=["N", "CA", "C", "O"], name=""):
    """Convert a PDB structure in to coordinates of target atoms from all AAs

    Args:
        struct: a Bio.PDB.Structure object representing the protein structure
        target_atoms: Target atoms which residues will be returned.
        name: String. Name of the structure

    Return:
        Dictionary with the pdb sequence, atom 3D coordinates and name.
    """
    output = {}
    # get AA sequence in the pdb structure
    pdb_seq = "".join(
        [three_to_one(res.get_resname()) for res in struct.get_residues()]
    )
    output["seq"] = pdb_seq
    # get the atom coords
    coords = np.asarray(
        [
            get_atom_coords(res, target_atoms=target_atoms)
            for res in struct.get_residues()
        ]
    )
    output["coords"] = coords.tolist()
    output["name"] = name
    return output


def parse_pdb_gz_to_json_record(parser, sequence, pdb_file_path, name=""):
    """
    Reads and reformats a pdb strcuture into a dictionary.

    Args:
        parser: a Bio.PDB.PDBParser or Bio.PDB.MMCIFParser instance.
        sequence: String. Sequence of the structure.
        pdb_file_path: String. Path to the pdb file.
        name: String. Name of the protein.

    Return:
        Dictionary with the pdb sequence, atom 3D coordinates and name.
    """
    struct = parse_pdb_structure(parser, sequence, pdb_file_path)
    record = structure_to_coords(struct, name=name)
    return record

if __name__=='__main__':
    args = parse_args()
    df = pd.read_csv(args.data_file)

    # PDB parser
    sloppyparser = PDBParser(
        QUIET=True,
        PERMISSIVE=True,
        structure_builder=xpdb.SloppyStructureBuilder(),
    )

    # 2. Parallel parsing structures and converting to protein records
    records = Parallel(n_jobs=-1)(
        delayed(parse_pdb_gz_to_json_record)(
            sloppyparser,
            df.iloc[i]["sequence"],
            df.iloc[i]["pdbpath"],
            df.iloc[i]["pdbpath"].split("/")[-1],
        )
        for i in tqdm(range(df.shape[0]))
    )
    
    target_col = args.target_col
    smiles_col = args.smiles_col
    
    for i, rec in enumerate(records):
        row = df.iloc[i]
        target = row[target_col]
        smiles = row[smiles_col]
        rec["target"] = target
        rec["ec"] = row["ec"]
        rec["taxonomy_id"] = row["taxonomy_id"]

    outprefix = args.out_prefix
    # 4. write to disk
    print("number of records:", len(records))
    outfile = os.path.join(f"{outprefix}_pdbrecords.json")
    json.dump(records, open(outfile, "w"))
    print('Success!')