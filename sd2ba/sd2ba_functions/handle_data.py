#!/usr/bin/env python3

from Bio.PDB import PDBParser

def get_solved_residues_from_pdb(pdb_code, pdb_file):
    """WIP"""

    structure = PDBParser().get_structure(pdb_code, pdb_file)

    # The Structure object follows the SMCRA (Structure/Model/Chain/Residue/Atom) architecture:
    model = structure[0]
    # Chain of the given reference PDB structure.
    chain = model[pdb_code[-1]]
    # Residues (Disordered residues not included) are classes that cotain: (id, resname, segid)
    residues = chain.get_residues()

    d3to1 = {
        "CYS": "C",
        "ASP": "D",
        "SER": "S",
        "GLN": "Q",
        "LYS": "K",
        "ILE": "I",
        "PRO": "P",
        "THR": "T",
        "PHE": "F",
        "ASN": "N",
        "GLY": "G",
        "HIS": "H",
        "LEU": "L",
        "ARG": "R",
        "TRP": "W",
        "ALA": "A",
        "VAL": "V",
        "GLU": "E",
        "TYR": "Y",
        "MET": "M",
        }

    solved_residues = {}

    for r in residues:
        # If the residue it's actually an aminoacid (could be water, an artefact...):
        if r.resname in d3to1:
            # "id" is a tuple that contains contains:
            # (hetero-flag[0], sequence identifier[1] and insertion code[2]).
            # "resname" is a 3 letter code of the name of the residue.
            # With d3to1 we convert the 3 letter to 1 letter code.
            solved_residues[r.id[1]] = {
                "resname": d3to1[r.resname],
                "insertion_code": r.id[2],
            }

    return solved_residues
