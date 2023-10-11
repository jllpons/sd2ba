#!/usr/bin/env python3

"""
sd2ba_fetch-data.py

From a list of UniProt IDs, this script fetches the AA (from UniProt) and
nucleotide (from ENA) sequences and writes them in two fasta files. FASTA
headers are in the format: >UniprotAccession_ENAAccession.

Accepts input from stdin or from a file with a list of UniProt IDs. Also
accepts a UniRef cluster ID as input.
"""

import argparse
import logging
import os
import sys

from sd2ba_functions.fetch_data import get_data_from_uniref, get_uniprot_entry_data, get_ena_nucleotide_sequence
from sd2ba_functions.output import write_fasta


class Protein:

    def __init__(self, input_id: str):
        self.input_id = input_id
        self.uniprot_accession = ""
        self.ena_accession = ""
        self.fasta_header = ""
        self.uniprot_aa_sequence = ""
        self.ena_nucleotide_sequence = ""


__version__ = "0.0.0"


def main():

    description = """
    Description: From a list of UniProt IDs, this script fetches the AA (from UniProt) and
    nucleotide (from ENA) sequences and writes them in two fasta files. FASTA
    headers are in the format: >UniprotAccession_ENAAccession.\n\n

    Accepts input from stdin or from a file with a list of UniProt IDs. Also
    accepts a UniRef cluster ID as input.
    """

    parser = argparse.ArgumentParser(
            description=description,
            usage="sd2ba_fetch-data.py [options]",
            )

    parser.add_argument(
            "-f",
            metavar="<file>",
            type=str,
            help="File with a list of Uniprot IDs separated by newlines. Mutually exclusive with --uniref",
            )
    parser.add_argument(
            "--uniref",
            metavar="<id>",
            type=str,
            help="UniRef cluster ID. Mutually exclusive with -f",
            )
    parser.add_argument(
            "-o", "--output",
            metavar="<path>",
            type=str,
            default=f"{os.getcwd()}/sd2ba_fetch-data_output",
            help="Output directory [Default: $CWD/sd2ba_fetch-data_output]",
            )

    parser.add_argument(
            "-V", "--version",
            action="version",
            version=f"%(prog)s {__version__}",
            )
    parser.add_argument(
            "-l", "--logging",
            metavar="<level>",
            type=str,
            default="INFO",
            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
            help="Set the logging level [Default: INFO] [Choices: DEBUG, INFO, WARNING, ERROR, CRITICAL]",
            )

    args = parser.parse_args()

    if args.uniref and args.f:
        parser.print_help()
        sys.exit("\n** ERROR: --uniref and -f are mutually exclusive **\n")

    if not args.uniref and not args.f:
        # If stdin is not connected to a terminal, (i.e. in a pipe) read from it
        if not sys.stdin.isatty():
            stdin = sys.stdin.read().splitlines()
            if len(stdin) == 0:
                parser.print_help()
                sys.exit("\n** ERROR: stdin is empty **\n")
        else:
            parser.print_help()
            sys.exit("\n** ERROR: runnig without --uniref or -f requires stdin **\n")

    output_dir_path = args.output_directory
    if not os.path.exists(output_dir_path):
        os.makedirs(output_dir_path)
    output_dir_path = os.path.abspath(output_dir_path)

    # Set up logger
    log_level = getattr(logging, args.logging.upper())
    logging.basicConfig(
                filename=f"{output_dir_path}/sd2ba_fetch-data.log",
                encoding="utf-8",
                format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                level=log_level,
                )

    logging.info(f"{__file__} {__version__}")
    logging.info("Script started")
    logging.info("CMD: " + " ".join(sys.argv))
    logging.debug(f"Output directory was created: {output_dir_path}")

    if args.uniref:
        # Accessing UniRef data and obtaining the list of proteins ids
        p_list = get_data_from_uniref(args.uniref)

    elif args.ids:
        # Reading the list of proteins ids
        with open(args.ids, "r") as f:
            p_list = f.read().splitlines()

    elif stdin:
        p_list = stdin

    # Creating a list of Protein objects from the list of proteins ids
    proteins = [Protein(input_id=i) for i in p_list]

    for p in proteins:
        uniprot_entry_data = get_uniprot_entry_data(p.input_id)

        if uniprot_entry_data["successful"] == True:
            p.uniprot_accession = uniprot_entry_data["uniprot_accession"]
            p.ena_accession = uniprot_entry_data["ena_accession"]
            p.uniprot_aa_sequence = uniprot_entry_data["uniprot_aa_sequence"]

            ena_nucleotide_sequence = get_ena_nucleotide_sequence(p.ena_accession)
            if ena_nucleotide_sequence["successful"] == True:
                p.ena_nucleotide_sequence = ena_nucleotide_sequence["fasta"]["sequence"]
                p.fasta_header = f">{p.uniprot_accession}_{p.ena_accession}"

            else:
                proteins.remove(p)
                logging.warning(
                        f"{p.input_id} have been removed "
                        + "from the analysis due to problems accesing its ENA entry. "
                        )

        else:
            proteins.remove(p)
            logging.warning(
                    f"{p.input_id} have been removed "
                    + "from the analysis due to problems accesing its Uniprot entry. "
                    )

    fasta_aa_sequences = []
    fasta_nucleotide_sequences = []
    for p in proteins:
        fasta_aa_sequences.append(p.fasta_header)
        fasta_aa_sequences.append(p.uniprot_aa_sequence)

        fasta_nucleotide_sequences.append(p.fasta_header)
        fasta_nucleotide_sequences.append(p.ena_nucleotide_sequence)

    if args.uniref:
        write_fasta(
                path=f"{output_dir_path}/{args.uniref}_aa_sequences.fasta",
                content="\n".join(fasta_aa_sequences),
                )
    else:
        write_fasta(
                path=f"{output_dir_path}/aa_sequences.fasta",
                content="\n".join(fasta_aa_sequences),
                )
    logging.info(f"The AA sequences from {len(proteins)} proteins present "
                 + f"in {args.uniref} were extracted from Uniprot and written to "
                 + f"{output_dir_path}/{args.uniref}_aa_sequences.fasta"
                 )

    if args.uniref:
        write_fasta(
                path=f"{output_dir_path}/{args.uniref}_nucleotide_sequences.fasta",
                content="\n".join(fasta_nucleotide_sequences),
                )
    else:
        write_fasta(
                path=f"{output_dir_path}/nucleotide_sequences.fasta",
                content="\n".join(fasta_nucleotide_sequences),
                )
    logging.info(f"The nucleotide sequences from {len(proteins)} proteins present "
                 + f"in {args.uniref} were extracted from ENA and written to "
                 + f"{output_dir_path}/{args.uniref}_nucleotide_sequences.fasta"
                 )

    logging.info("Script finished successfully")


if __name__ == "__main__":
    main()
