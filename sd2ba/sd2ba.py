#!/usr/bin/env python3

"""
sd2ba: Segmented Domain to Breakpoint Analysis
"""

# TODO: Remove later
import pdb

import argparse
from functools import total_ordering
import logging
import os
import sys

from requests.adapters import proxy_from_url

from sd2ba_functions.fetch_data import get_data_from_uniref, get_uniprot_entry_data, get_pdb_file
from sd2ba_functions.handle_data import get_solved_residues_from_pdb
from sd2ba_functions.output import write_fasta


__version__ = "0.0"


@total_ordering
class Protein:

    def __init__(self, input_id: str):
        self.input_id = input_id
        self.uniprot_accession = ""
        self.ena_accession = ""
        self.uniprot_aa_sequence = ""
        self.similarity_value = 0

    def __eq__(self, other):
        return self.similarity_value == other.similarity_value

    def __lt__(self, other):
        return self.similarity_value < other.similarity_value

    def __gt__(self, other):
        return self.similarity_value > other.similarity_value


class ReferenceProtein(Protein):

    def __init__(self, input_id: str, aa_with_known_positions: dict):
        super().__init__(input_id)
        self.aa_with_known_positions = aa_with_known_positions


def main():

    description = """sd2ba: Segmented Domain to Breakpoint Analysis"""

    parser = argparse.ArgumentParser(
                        description=description,
                        usage="sd2ba.py PDB-CODE UNIREF-CODE [options]"
                        )

    parser.add_argument(
            "pdb_code",
            metavar="PDB-CODE",
            type=str,
            help="PDB code of the reference protein including the chain",
            )
    parser.add_argument(
            "uniref_code",
            metavar="UNIREF-CODE",
            type=str,
            help="TESTING"
            )

    parser_output = parser.add_argument_group("output options")
    parser_output.add_argument(
                    "-o", "--output-directory",
                    metavar="STR",
                    type=str,
                    default="sd2ba_output",
                    help = "path/to/output_directory [Default: s2ba_output]"
                    )

    parser.add_argument(
                    "-l", "--logging",
                    metavar="STR",
                    type=str,
                    default="INFO",
                    help = "Specify the logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL. [Default: INFO]"
                    )
    parser.add_argument(
                    "-v", "--version",
                    action="version",
                    version=f"sd2ba.py v{__version__}"
            )

    args = parser.parse_args()

    # PDB code
    pdb_code = args.pdb_code.upper()
    if len(pdb_code) != 5:
        sys.exit(f"\n** The PDB code of the reference protein must specify the chain and provided was: {pdb_code} **")

    # Output directory
    output_path = args.output_directory
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = os.path.abspath(output_path)

    # Set up logger
    log_level = getattr(logging, args.logging.upper())
    logging.basicConfig(
                filename=f"{output_path}/sd2ba.log",
                encoding="utf-8",
                format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                level=log_level,
                )

    logging.info(f"sd2ba.py v{__version__}")
    logging.info("CMD: " + " ".join(sys.argv))
    logging.debug(f"Output directory was created: {output_path}")

    uniref_data = get_data_from_uniref(args.uniref_code)
    proteins = []
    for i in uniref_data:
        proteins.append(Protein(input_id=i))

    for p in proteins:
        uniprot_entry_data = get_uniprot_entry_data(p.input_id)

        if uniprot_entry_data["successful"] == True:
            p.uniprot_accession = uniprot_entry_data["uniprot_accession"]
            p.ena_accession = uniprot_entry_data["ena_accession"]
            p.uniprot_aa_sequence = uniprot_entry_data["uniprot_aa_sequence"]

        else:
            proteins.remove(p)
            logging.warning(
                    f"{p.input_id} present in {args.uniref_code} have been removed "
                    + "from the analysis due to problems accesing its uniprot entry. "
                    )

    fasta_uniprot_aa_sequences = []
    for p in proteins:
        fasta_uniprot_aa_sequences.append(f">{p.uniprot_accession}_{p.ena_accession}")
        fasta_uniprot_aa_sequences.append(p.uniprot_aa_sequence)

    write_fasta(f"{output_path}/uniprot_aa_sequences.fasta", "\n".join(fasta_uniprot_aa_sequences))

    # Requesting PDB file to RSCB-PDB.
    # Any unsuccessful request will result in critical error and the exit of the script.
    # The PDB file is saved in the output directory
    pdb_file = output_path + f"/{pdb_code}.pdb"
    with open(pdb_file, "w") as handle:
        handle.write(get_pdb_file(pdb_code[:4]))
    logging.debug(f"PDB file for {pdb_code} was saved in {pdb_file}")

    ref_protein = ReferenceProtein(
            input_id=pdb_code,
            aa_with_known_positions=get_solved_residues_from_pdb(pdb_code, pdb_file)
            )

    for p in proteins:
            p.similarity_value = calculate_similarity(ref_protein, p)


if __name__ == "__main__":
    main()
