#!/usr/bin/env python3

"""
Author:         Joan Lluis Pons Ramon

Usage:          python3 id_aant.py [options]

Description:    From a list of UniProt IDs, this script fetches the amino acid
                (from UniProt) and nucleotide (from ENA) sequence and writes
                them in two fasta files.
                FASTA headers will appear as: >UniprotAccession_ENAAccession.
                Accepts IDs list from sdin or from file if -f is used.
"""

import json
import logging
import os
import sys

import requests
from requests.adapters import HTTPAdapter, Retry

from sd2ba_functions.API_URL import UNIPROT_ENTRY_JSON_API_URL, ENA_NUCLEOTIDE_SEQUENCE_API_URL


__version__ = "0.0.0"


class Arguments:
    """
    Arguments object with the parsed arguments. Default values are set in the
    __init__ method.
    """

    def __init__(self):
        self.id_file = None
        self.output_directory = f"{os.getcwd()}/id_aant_output"
        self.log_level = "INFO"


class Cli:
    """
    Command line interface object. It parses the arguments provided by the user
    and prints the usage, description and help messages if needed.
    """

    def __init__(self):

        self.usage = "Usage: python3 id_aant.py [options]"

        self.description = """
        Description: From a list of UniProt IDs, this script fetches the
        amino acid (from UniProt) and nucleotide (from ENA) sequence and writes
        them in two fasta files.
        FASTA headers will appear as: >UniprotAccession_ENAAccession.
        Accepts IDs list from sdin or from file if -f is used.
        """

        self.help_msg = """
        -f    <file>    File with a list of UniProt IDs separated by newlines.
        -o    <path>    Output directory [Default: $CWD/id_aant_output]
        -l    <level>   Set the logging level [Default: INFO] [Choices: DEBUG, INFO, WARNING, ERROR, CRITICAL]
        -V              Show version and exit
        -h, --help      Show this help message and exit
        """

    def print_help_and_exit(self) -> None:
        """
        Prints the help message and exits
        """
        print(self.usage)
        print(self.description)
        print(self.help_msg)
        sys.exit(1)


    def chek_filepath(self, path: str) -> bool:
        """
        Checks if the path is valid and returns the absolute file_path

        Args:
            path (str): The path to check.

        Returns:
            bool: True if the path is valid.
        """
        if not os.path.exists(path):
            sys.exit(f"\nError: {path} does not exist\n")
        if not os.path.isfile(path):
            sys.exit(f"\nError: {path} is not a file\n")
        return True


    def parse_arguments(self, arguments: list) -> object:
        """
        Parses the arguments provided by the user.

        Args:
            arguments (list): The arguments provided by the user.

        Returns:
            Arguments (object): An Arguments object with the parsed arguments.
        """
        if "-h" in arguments or "--help" in arguments:
            self.print_help_and_exit()

        parsed_args = Arguments()
        for indx, arg in enumerate(arguments):
            if arg == "-f":
                try:
                    file_path = arguments[indx + 1]
                except IndexError:
                    sys.exit("\nError: -f requires a file path\n")
                if self.chek_filepath(file_path):
                   parsed_args .id_file = file_path

            elif arg == "-o":
                try:
                    output_directory = arguments[indx + 1]
                except IndexError:
                    sys.exit("\nError: -o requires a path\n")
                parsed_args.output_directory = output_directory

            elif arg == "-l":
                try:
                    log_level = arguments[indx + 1].upper()
                except IndexError:
                    sys.exit("\nError: -l requires a log level\n")
                if log_level not in ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]:
                    sys.exit("\nError: log level must be one of DEBUG, INFO, WARNING, ERROR, CRITICAL\n")
                parsed_args.log_level = log_level

            elif arg == "-V":
                print(f"{__file__} {__version__}")
                sys.exit(0)

        return parsed_args


class Protein:
    """
    Protein object with the data fetched from UniProt and ENA.
    """

    def __init__(self, input_id: str):
        self.input_id = input_id
        self.uniprot_accession = ""
        self.ena_accession = ""
        self.uniprot_aa_sequence = ""
        self.ena_nucleotide_sequence = ""

    def generate_fasta_header(self) -> str:
        """
        Generates the fasta header for the protein.

        Returns:
            str: The fasta header.
        """
        return f">{self.uniprot_accession}_{self.ena_accession}"


class FetchDataResult:
    """
    FetchDataResult object with the result of the request.
    Can be successful or unsuccessful.
    """

    def __init__(self, successful: bool, data: dict):
        self.successful = successful
        self.data = data


class FetchData:
    """
    FetchData object with the methods to fetch data from UniProt and ENA.
    """

    def __init__(self):
        pass


    def make_request_with_retries(self, url: str, retries=3, delay=1) -> requests.Response:
        """
        Make a request to a given URL with retries in case of failure.

        Parameters:
            url (str): The URL to send the request to.
            retries (int): The maximum number of retries in case of failure (default: 3).
            delay (float): The delay in seconds between retries. Increses exponentially.(default: 1).

        Returns:
            requests.Response: The response object.
        """

        retry_strategy = Retry(
                            total=retries,
                            status_forcelist=[429,500,502,503,504],
                            allowed_methods=["GET", "POST"],
                            backoff_factor=delay,
                            )

        adapter = HTTPAdapter(max_retries=retry_strategy)
        http = requests.Session()
        http.mount("https://", adapter)

        response = http.get(url)
        return response


    def handle_uniprot_entry_json(self, j: str) -> FetchDataResult:
        """
        Handles the JSON response from a UniProt entry.

        Args:
            json (str): The JSON response from a UniProt entry.

        Returns:
            FetchDataResult: A FetchDataResult object with the result of the request.
                             FetchDataResult can be successful or unsuccessful.
        """

        response = json.loads(j)

        try:

            if response["entryType"] == "Inactive":
                return FetchDataResult(successful=False, data={})

            # UniProt accession code.
            uniprot_accession = response["primaryAccession"]

            for i in range(len(response["uniProtKBCrossReferences"])):

                if response["uniProtKBCrossReferences"][i]["database"] == "EMBL":

                    embl_data = response["uniProtKBCrossReferences"][i]

                    for j in range(len(embl_data["properties"])):

                        if embl_data["properties"][j]["key"] == "ProteinId":
                            # ENA accession code.
                            ena_accession = embl_data["properties"][j]["value"].split(".")[0]

            # Amino acid sequence.
            uniprot_aa_sequence = response["sequence"]["value"]

            # FIXME: some entries present errors
            if len(uniprot_accession) and len(ena_accession) and len(uniprot_aa_sequence):

                return FetchDataResult(
                        successful=True,
                        data={
                            "uniprot_accession": uniprot_accession,
                            "ena_accession": ena_accession,
                            "uniprot_aa_sequence": uniprot_aa_sequence,
                            }
                        )

            else:
                return FetchDataResult(successful=False, data={})

        except:
            return FetchDataResult(successful=False, data={})


    def get_uniprot_entry_data(self, code: str) -> FetchDataResult:
        """
        Fetches the data from a UniProt entry.

        Args:
            code (str): The UniProt accession code.

        Returns:
            FetchDataResult: A FetchDataResult object with the result of the request.
                             FetchDataResult can be successful or unsuccessful.
        """

        url = UNIPROT_ENTRY_JSON_API_URL.format(code=code)

        response = self.make_request_with_retries(url)

        if response.status_code == 200:
            logging.debug(f"UniProt request for {code} entry was successful. "
                          + f"Url used was {url}.")

            try:
                entry_data = self.handle_uniprot_entry_json(response.text)

                if entry_data.successful == True:

                    return FetchDataResult(successful=True, data=entry_data.data)

                logging.info(
                    f"Error in handling UniProt JSON response for {code}. "
                    + "The structure of the JSON response may have changed. "
                    + "Check FetchData.handle_uniprot_entry_json function or "
                    + "run 'python -m unittest test.test_id_aant'."
                    )

                return FetchDataResult(successful=False, data={})

            except:
                logging.info(
                    f"Error in handling UniProt JSON response for {code}. "
                    + "The structure of the JSON response may have changed. "
                    + "Check FetchData.handle_uniprot_entry_json function or "
                    + "run 'python -m unittest test.test_id_aant'."
                        )
                return FetchDataResult(successful=False, data={})

        logging.info(f"UniProt request for {code} entry was unsuccessful. "
                     + f"Url used was {url}. Response status code  "
                     + f"was {response.status_code}.")
        return FetchDataResult(successful=False, data={})


    def get_ena_nucleotide_sequence(self, code: str) -> FetchDataResult:
        """
        Fetches the nucleotide sequence from ENA for a given accession code.

        Args:
            code (str): The ENA accession code.

        Returns:
            FetchDataResult: A FetchDataResult object with the result of the request.
                             FetchDataResult can be successful or unsuccessful.
        """

        url = ENA_NUCLEOTIDE_SEQUENCE_API_URL.format(code=code)

        response = self.make_request_with_retries(url)

        if response.status_code == 200:
            logging.debug(
                    f"ENA request for {code} nucleotide sequence was successful. "
                    + f"Url used was {url}.")

            ena_fasta = response.text

            if ena_fasta.startswith(">"):
                ena_fasta = ena_fasta.split("\n")
                ena_fasta = {
                        "header": ena_fasta[0].strip(),
                        "sequence": "".join(ena_fasta[1:]).strip("\n"),
                        }

            else:
                logging.info(
                        f"Fast sequence for {code} nucleotide sequence was not found "
                        + "or was not recognized as a fasta sequence. "
                        + f"Url used was {url}."
                        )
                return FetchDataResult(successful=False, data={})

            return FetchDataResult(successful=True, data=ena_fasta)

        logging.info(
                f"ENA request for {code} nucleotide sequence was unsuccessful. "
                + f"Url used was {url}. "
                + f"Response status code was {response.status_code}.")
        return FetchDataResult(successful=False, data={})


def main():

    args = Cli().parse_arguments(sys.argv[1:])

    # IDs
    if args.id_file:
        with open(args.id_file, "r") as f:
            id_list = f.read().splitlines()
        if len(id_list) == 0:
            sys.exit("\nError: {args.id_file} is empty\n")
    else:
        if not sys.stdin.isatty():
            id_list = sys.stdin.read().splitlines()
            if len(id_list) == 0:
                sys.exit("\nError: no IDs were provided\n")
        else:
            sys.exit("\nError: runnig the script without -f requires IDs from stdin\n")

    # Output directory
    output_dir_path = args.output_directory
    if not os.path.exists(output_dir_path):
        os.makedirs(output_dir_path)
    output_dir_path = os.path.abspath(output_dir_path)

    # Logging
    logging.basicConfig(
                filename=f"{output_dir_path}/id_aant.log",
                encoding="utf-8",
                format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                level=args.log_level,
                )
    logging.info(f"{__file__} v{__version__}")
    logging.info("Script started")
    logging.info("CMD: " + " ".join(sys.argv))
    logging.debug(f"Output directory was created: {output_dir_path}")
    logging.info(f"{len(id_list)} IDs were provided")

    # Data fetching
    proteins = [Protein(input_id=i) for i in id_list]
    for p in proteins:
        uniprot_entry_data = FetchData().get_uniprot_entry_data(p.input_id)

        if uniprot_entry_data.successful == True:
            p.uniprot_accession = uniprot_entry_data.data["uniprot_accession"]
            p.ena_accession = uniprot_entry_data.data["ena_accession"]
            p.uniprot_aa_sequence = uniprot_entry_data.data["uniprot_aa_sequence"]

            ena_nucleotide_sequence = FetchData().get_ena_nucleotide_sequence(p.ena_accession)
            if ena_nucleotide_sequence.successful == True:
                p.ena_nucleotide_sequence = ena_nucleotide_sequence.data["sequence"]

            else:
                proteins.remove(p)
                logging.warning(f"{p.input_id} have been removed "
                                + "from the analysis due to problems accesing "
                                + "its Uniprot entry. ")

        else:
            proteins.remove(p)
            logging.warning(f"{p.input_id} have been removed "
                            + "from the analysis due to problems accesing "
                            + "its Uniprot entry. ")

    logging.info(f"{len(proteins)} proteins were fetched successfully")

    # Writing fasta files
    if len(proteins) > 0:
        fasta_aa_sequences = []
        fasta_nucleotide_sequences = []
        for p in proteins:
            fasta_aa_sequences.append(p.generate_fasta_header())
            fasta_aa_sequences.append(p.uniprot_aa_sequence)

            fasta_nucleotide_sequences.append(p.generate_fasta_header())
            fasta_nucleotide_sequences.append(p.ena_nucleotide_sequence)

        with open(f"{output_dir_path}/aa_sequences.fasta", "w") as f:
            f.write("\n".join(fasta_aa_sequences))
        with open(f"{output_dir_path}/nucleotide_sequences.fasta", "w") as f:
            f.write("\n".join(fasta_nucleotide_sequences))

        logging.info("Script finished successfully")
        sys.exit(0)

    sys.exit("\nError: no proteins were fetched successfully\n")


if __name__ == "__main__":
    main()
