#!/usr/bin/env python3

# TODO: Remove later
import pdb

import json
import logging

from Bio.PDB import PDBParser


logger = logging.getLogger(__name__)


def handle_uniref_json(response):
    """
    Handles the JSON response from UniRef and returns a list of UniRef IDs
    for the members of the cluster.

    Args:
        response (str): The JSON response from UniRef.

    Returns:
        list: A list of Uniprot IDs for the members of the cluster."
    """

    response = json.loads(response)

    members = response["members"]

    uniref_data = []
    for i in members:

        if i["memberIdType"] == "UniProtKB ID":
            uniref_data.append(i["accessions"][0])

        # TODO: if uniref it's the final option. Handle uniparc accessions.
        # elif i[] == Uniparc...

    logging.info(
            f"Number of members found in the UniRef cluster was {len(members)}. "
            + f"Number of proteins IDs included in the analysis was {len(uniref_data)}."
                  )

    return uniref_data

def handle_uniprot_entry_json(response):
    """
    Handles the JSON response from UniProt for an entry.

    Args:
        response (str): The JSON response string from UniProt.

    Returns:
        dict: A dictionary containing the extracted entry data. The dictionary has the following keys:
            - successful (bool): Indicates whether the handling was successful.
            - uniprot_accession (str): The UniProt accession code for the entry.
            - ena_accession (str): The ENA (European Nucleotide Archive) accession code for the entry.
            - uniprot_aa_sequence (str): The amino acid sequence of the entry.

    Notes:
        This function takes the JSON response from UniProt and extracts the relevant information.
        It retrieves the primary UniProt accession code, the ENA accession code from the EMBL database,
        and the amino acid sequence of the entry.
        If all the necessary data is successfully extracted, it is returned in a dictionary format with "successful" set to True.
        If any errors occur during the extraction or if the required data is missing, the function returns a dictionary with "successful" set to False.

    Raises:
        None.
    """

    response = json.loads(response)

    try:
        uniprot_accession = response["primaryAccession"]

        for i in range(len(response["uniProtKBCrossReferences"])):

            if response["uniProtKBCrossReferences"][i]["database"] == "EMBL":

                embl_data = response["uniProtKBCrossReferences"][i]

                for j in range(len(embl_data["properties"])):

                    if embl_data["properties"][j]["key"] == "ProteinId":
                        ena_accession = embl_data["properties"][j]["value"].split(".")[0]

        uniprot_aa_sequence = response["sequence"]["value"]

        # FIXME: some entries present errors
        if len(uniprot_accession) and len(ena_accession) and len(uniprot_aa_sequence):

            return {
                    "successful": True,
                    "uniprot_accession": uniprot_accession,
                    "ena_accession": ena_accession,
                    "uniprot_aa_sequence": uniprot_aa_sequence,
                    }

        else:
            return {"successful": False}

    except:
        return {"successful": False}


def get_solved_residues_from_pdb(pdb_code, pdb_file):
    """
    Retrieves the solved residues from a PDB file.

    Args:
        pdb_code (str): The PDB code of the structure.
        pdb_file (str): The path to the PDB file.

    Returns:
        dict: A dictionary containing information about the solved residues.
            The dictionary has the residue sequence identifier as the key and a
            dictionary with the following keys as the value:
                - 'resname': The one-letter code of the residue.
                - 'insertion_code': The insertion code of the residue.
    """


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
