#!/usr/bin/env python3


"""
sd2ba: Segmented Domain to Breakpoint Analysis
"""


import argparse
from functools import total_ordering
import logging
import os
import sys

from sd2ba_functions.fetch_data import (
        get_ena_nucleotide_sequence, get_pdb_file, get_hmm_file,
        get_aa_sequence_from_pdb, get_up_code_from_pdb_code, get_uniprot_entry_data
        )
from sd2ba_functions.handle_data import get_solved_residues_from_pdb, read_multiple_fasta
from sd2ba_functions.SCRIPT_ARGS import (
        HMMER_ARGS,
        SED_ARGS_AFA_TO_FASTA
        )
from sd2ba_functions.subprocess import run_subprocess


__version__ = "0.0.0"


@total_ordering
class Protein:

    def __init__(self, header: str):
        self.header = header
        self.aa_sequence = ""
        self.nt_sequence = ""
        self.aligned_aa_sequence = ""
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
        self.aa_sequence_pdb = ""
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
    try:
        # get_solved_residues_from_pdb() accesses the PDB file
        # and returns a dictionary with the solved residues
        ref_protein.aa_with_known_positions=get_solved_residues_from_pdb(pdb_code, pdb_file)
    except:
        logging.critical(f"Error while parsing the PDB file of {pdb_code}")
        logging.critical("script ended")
        sys.exit("\n** Error while parsing the PDB file, check the log file for more information **")

    ref_protein.aa_sequence_pdb = get_aa_sequence_from_pdb(pdb_code[:4])
    uniprot_entry = get_uniprot_entry_data(get_up_code_from_pdb_code(pdb_code[:4]))
    if uniprot_entry["successful"] == False:
        logging.critical(f"Error accessing the Uniprot entry of {pdb_code}")
        logging.critical("script ended")
        sys.exit("\n** Error accessing the Uniprot entry, check the log file for more information **")
    ref_protein.aa_sequence = uniprot_entry["uniprot_aa_sequence"]

    ena_entry = get_ena_nucleotide_sequence(uniprot_entry["ena_accession"])
    if ena_entry["successful"] == False:
        logging.critical(f"Error accessing the ENA entry of {pdb_code}")
        logging.critical("script ended")
        sys.exit("\n** Error accessing the ENA entry, check the log file for more information **")
    ref_protein.nt_sequence = ena_entry["fasta"]["sequence"]

    logging.info(f"Reference protein {pdb_code} has {len(ref_protein.aa_with_known_positions)} solved residues: "
            + f"[{list(ref_protein.aa_with_known_positions.keys())[0]}"
            + f"..{list(ref_protein.aa_with_known_positions.keys())[-1]}]"
            )

    hmm_file = output_path + f"/{pfam_code}.hmm"
    with open(hmm_file, "w") as handle:
        # get_hmm_file() returns a string with the obtained HMM file
        # If the request was unsuccessful, the script will exit
        handle.write(get_hmm_file(pfam_code))
    logging.debug(f"HMM file for {pfam_code} was saved in {hmm_file}")

    # Parsing the amino and nucleotide fasta files
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

    # The number of amino acid and nucleotide sequences must be the same
    if len(aa_sequences) != len(nt_sequences):
        logging.critical("The number of amino acid and nucleotide sequences must be the same")
        logging.critical("script ended")
        sys.exit("\n** The number of amino acid and nucleotide sequences must be the same **")

    # Loading the amino acid and nucleotide sequences into Protein objects
    proteins = []
    for header in aa_sequences.keys():
        p = Protein(header=header)
        p.aa_sequence = aa_sequences[header]
        p.nt_sequence = nt_sequences[header]
        proteins.append(p)

    # The reference protein will be the first one in the list
    proteins.insert(0, ref_protein)

    input_with_pdb_aa_fasta = output_path + f"/input_with_{pdb_code}_amino.fasta"
    input_with_pdb_nt_fasta = output_path + f"/input_with_{pdb_code}_nucleotide.fasta"

    # Rerwiting the input files with the reference protein as the first one
    with open(input_with_pdb_aa_fasta, "w") as handle:
        for p in proteins:
            handle.write(f">{p.header}\n{p.aa_sequence}\n")

    with open(input_with_pdb_nt_fasta, "w") as handle:
        for p in proteins:
            handle.write(f">{p.header}\n{p.nt_sequence}\n")

    logging.info(
            f"Nuclotide and amino acid sequences file were rewritten with {pdb_code} as the first sequence "
            + f"and saved in {input_with_pdb_aa_fasta} and {input_with_pdb_nt_fasta}"
            )

    # Running HMMER to make a multiple sequence alignment
    hmmalign_stdout = run_subprocess(
                        "HMMALIGN",
                        HMMER_ARGS
                            .format(
                                hmm_filepath=hmm_file,
                                aa_seq_filepath=input_with_pdb_aa_fasta,
                                )
                            .split(),
                            )
    logging.info("HMMER's HMMALIGN  was used to make a multiple sequence alignment")

    # The hmmalign output is written in afa format. We will convert it to fasta
    # and save it to a file
    hmmalign_output = output_path + "/hmmalign_output.fasta"
    with open(hmmalign_output, "w") as handle:
        sed_stdout = run_subprocess(
                        "SED",
                        SED_ARGS_AFA_TO_FASTA,
                        input=hmmalign_stdout,
                        )
        handle.write(sed_stdout.decode("utf-8"))

    logging.info("SED was called to convert the hmmalign output to fasta format")
    logging.info(f"HMMER output was saved in {hmmalign_output}")

    # Loading the multiple sequence alignment sequences into the Protein objects
    msa_proteins = []
    with open(hmmalign_output, "r") as handle:
        records = read_multiple_fasta(handle.read())

        if records["successful"] == False:
            logging.critical(f"Error while parsing the multiple sequence alignment file {hmmalign_output}")
            logging.critical("script ended")
            sys.exit("\n** Error while parsing the multiple sequence alignment file, check the log file for more information **")

        msa_proteins = records["data"]

    logging.info(f"{len(msa_proteins)} sequences were read from {hmmalign_output}")

    for p in proteins:
        p.aligned_aa_sequence = msa_proteins[p.header]


if __name__ == "__main__":
    main()
