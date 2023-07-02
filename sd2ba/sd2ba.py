#!/usr/bin/env python3

"""
sd2ba: Segmented Domain to Breakpoint Analysis
"""

import argparse
from functools import total_ordering
import logging
import os
import sys

from sd2ba_functions.fetch_data import get_pdb_file
from sd2ba_functions.handle_data import get_solved_residues_from_pdb


__version__ = "0.0"


@total_ordering
class Protein:

    def __init__(self, input_id: str):
        self.input_id = input_id
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


    # Requesting PDB file to RSCB-PDB.
    # Any unsuccessful request will result in critical error and the exit of the script.
    # The PDB file is saved in the output directory
    pdb_file = output_path + f"/{pdb_code}.pdb"
    with open(pdb_file, "w") as handle:
        handle.write(get_pdb_file(pdb_code[:4]).text)
    logging.debug(f"PDB file for {pdb_code} was saved in {pdb_file}")

    ref_protein = ReferenceProtein(
            input_id=pdb_code,
            aa_with_known_positions=get_solved_residues_from_pdb(pdb_code, pdb_file)
            )


if __name__ == "__main__":
    main()
