#!/usr/bin/env python3


"""
sd2ba: Segmented Domain to Breakpoint Analysis
"""


import argparse
from functools import total_ordering
import logging
import os
import sys

from sd2ba_functions.fetch_data import get_pdb_file, get_hmm_file
from sd2ba_functions.handle_data import get_solved_residues_from_pdb, read_multiple_fasta


__version__ = "0.0.0"


@total_ordering
class Protein:

    def __init__(self, header: str):
        self.header = header
        self.aa_sequence = ""
        self.nt_sequence = ""
        self.similarity_value = 0

    def __eq__(self, other):
        return self.similarity_value == other.similarity_value

    def __lt__(self, other):
        return self.similarity_value < other.similarity_value

    def __gt__(self, other):
        return self.similarity_value > other.similarity_value


class ReferenceProtein(Protein):

    def __init__(self, header: str):
        super().__init__(header)
        self.aa_with_known_positions = {}


def main():

    description = """sd2ba.py is a tool to identify the breakpoint of a segmented domain in a protein family"""

    parser = argparse.ArgumentParser(
                        description=description,
                        usage="sd2ba.py PDB PFAM amino.fasta nucleotide.fasta [options]"
                        )

    parser.add_argument(
            "pdb_code",
            metavar="PDB",
            type=str,
            help="PDB code of the reference protein including the chain",
            )
    parser.add_argument(
            "pfam_code",
            metavar="PFAM",
            type=str,
            help="PFAM code of the domain of interest. The associated HMM will be downloaded.",
            )
    parser.add_argument(
            "aa_fasta",
            metavar="amino.fasta",
            type=str,
            help="Fasta file with the amino acid sequences of the proteins of interest.",
            )
    parser.add_argument(
            "nt_fasta",
            metavar="nucleotide.fasta",
            type=str,
            help="Fasta file with the nucleotide sequences of the proteins of interest.",
            )

    parser_output = parser.add_argument_group("output options")
    parser_output.add_argument(
                    "-o", "--output-directory",
                    metavar="STR",
                    type=str,
                    default="sd2ba_output",
                    help = "Output directory [Default: $CWD/s2ba_output]"
                    )

    parser.add_argument(
                    "-l", "--logging",
                    metavar="STR",
                    type=str,
                    default="INFO",
                    help = "Set the logging level [Default: INFO] [Choices: DEBUG, INFO, WARNING, ERROR, CRITICAL]"
                    )
    parser.add_argument(
                    "-V", "--version",
                    action="version",
                    version=f"%(prog)s v{__version__}",
            )

    args = parser.parse_args()

    # PDB code
    pdb_code = args.pdb_code.upper()
    if len(pdb_code) != 5:
        sys.exit(f"\n** The PDB code of the reference protein must specify the chain and provided was: {pdb_code} **")
    # PFAM code
    pfam_code = args.pfam_code.upper()

    # Amino fasta file
    aa_fasta = args.aa_fasta
    if not os.path.isfile(aa_fasta):
        sys.exit(f"\n** The provided amino fasta file does not exist: {aa_fasta} **")
    # Nucleotide fasta file
    nt_fasta = args.nt_fasta
    if not os.path.isfile(nt_fasta):
        sys.exit(f"\n** The provided nucleotide fasta file does not exist: {nt_fasta} **")

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

    logging.info(f"{__file__} v{__version__}")
    logging.info("script started")
    logging.info("CMD: " + " ".join(sys.argv))
    logging.debug(f"Output directory was created: {output_path}")

    # Requesting PDB file to RSCB-PDB.
    # Any unsuccessful request will result in critical error and the exit of the script.
    # The PDB file is saved in the output directory
    pdb_file = output_path + f"/{pdb_code}.pdb"
    with open(pdb_file, "w") as handle:
        # get_pdb_file() returns a string with the obtained PDB file
        # If the request was unsuccessful, the script will exit
        handle.write(get_pdb_file(pdb_code[:4]))
    logging.debug(f"PDB file for {pdb_code} was saved in {pdb_file}")

    ref_protein = ReferenceProtein(header=pdb_code)
    # get_solved_residues_from_pdb() accesses the PDB file
    # and returns a dictionary with the solved residues
    try:
        ref_protein.aa_with_known_positions=get_solved_residues_from_pdb(pdb_code, pdb_file)
    except:
        logging.critical(f"Error while parsing the PDB file of {pdb_code}")
        logging.critical("script ended")
        sys.exit("\n** Error while parsing the PDB file, check the log file for more information **")

    logging.info(f"Reference protein {pdb_code} has {len(ref_protein.aa_with_known_positions)} solved residues")

    hmm_file = output_path + f"/{pfam_code}.hmm"
    with open(hmm_file, "w") as handle:
        # get_hmm_file() returns a string with the obtained HMM file
        # If the request was unsuccessful, the script will exit
        handle.write(get_hmm_file(pfam_code))
    logging.debug(f"HMM file for {pfam_code} was saved in {hmm_file}")

    aa_sequences = {}
    with open(aa_fasta, "r") as handle:
        records = read_multiple_fasta(handle.read())

        if records["successful"] == False:
            logging.critical(f"Error while parsing the amino fasta file {aa_fasta}")
            logging.critical("script ended")
            sys.exit("\n** Error while parsing the amino fasta file, check the log file for more information **")

        aa_sequences = records["data"]

    logging.info(f"{len(aa_sequences)} amino acid sequences were read from {aa_fasta}")

    nt_sequences = {}
    with open(nt_fasta, "r") as handle:
        records = read_multiple_fasta(handle.read())

        if records["successful"] == False:
            logging.critical(f"Error while parsing the nucleotide fasta file {nt_fasta}")
            logging.critical("script ended")
            sys.exit("\n** Error while parsing the nucleotide fasta file, check the log file for more information **")

        nt_sequences = records["data"]

    logging.info(f"{len(nt_sequences)} nucleotide sequences were read from {nt_fasta}")

    if len(aa_sequences) != len(nt_sequences):
        logging.critical("The number of amino acid and nucleotide sequences must be the same")
        logging.critical("script ended")
        sys.exit("\n** The number of amino acid and nucleotide sequences must be the same **")

    proteins = []
    for header in aa_sequences.keys():
        p = Protein(header=header)
        p.aa_sequence = aa_sequences[header]
        p.nt_sequence = nt_sequences[header]
        proteins.append(p)


if __name__ == "__main__":
    main()
