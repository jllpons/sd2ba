#!/usr/bin/env python3

import logging
import sys

import requests

from requests.adapters import HTTPAdapter, Retry

from sd2ba_functions.api_urls import UNIPROT_UNIREF_JSON_API_URL, RSCB_PDB_API_URL, UNIPROT_ENTRY_JSON_API_URL
from sd2ba_functions.handle_data import handle_uniref_json, handle_uniprot_entry_json


logger = logging.getLogger(__name__)


def make_request_with_retries(url, retries=3, delay=1):
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


def get_data_from_uniref(code):
    """
    Accesses uniprot and obtains all of the uniprot accession codes from a given
    uniref code
    """

    url = UNIPROT_UNIREF_JSON_API_URL.format(code=code)

    response = make_request_with_retries(url)

    if response.status_code == 200:
        logging.debug(
                f"UniProt-UniRef request for {code} was successful. "
                + f"Url used was {url}."
                )
        try:
            return handle_uniref_json(response.text)
        except:
            logging.critical(
                    f"Error in handling UniProt JSON response for {code}. "
                    + "The structure of the JSON response may have changed. "
                    + "Check handle_uniref_json function in sd2ba_functions.handle_data"
                    )
            sys.exit("\n** CRITICAL: Error in handling UniProt response. Check log file for detailed info. **")
    else:
        logging.critical(
                f"UniProt request for {code} was unsuccessful. "
                + f"Url used was {url}. The script has stopped."
                )
        sys.exit("\n** CRITICAL: Request to UniProt failed. Check log file for detailed info. **")


def get_uniprot_entry_data(code):
    """
    Retrieves UniProt entry data for a given code.

    Args:
        code (str): The UniProt code for the entry.

    Returns:
        dict: A dictionary containing the retrieved entry data. The dictionary has the following keys:
            - successful (bool): Indicates whether the retrieval was successful.
            - uniprot_accession (str): The UniProt accession code for the entry.
            - ena_accession (str): The ENA (European Nucleotide Archive) accession code for the entry.
            - uniprot_aa_sequence (str): The amino acid sequence of the entry.

    Notes:
        This function makes a request to the UniProt API using the provided code and retrieves the corresponding entry data.
        If the retrieval is successful, the relevant information is extracted and returned in a dictionary format.
        If any errors occur during the retrieval or handling of the response, the function logs appropriate messages and returns a dictionary with "successful" set to False.

    Raises:
        None.
    """

    url = UNIPROT_ENTRY_JSON_API_URL.format(code=code)

    response = make_request_with_retries(url)

    if response.status_code == 200:
        logging.debug(
                f"UniProt request for {code} entry was successful. "
                + f"Url used was {url}."
                )

        try:
            entry_data = handle_uniprot_entry_json(response.text)

            if entry_data["successful"] == True:
                return {
                        "successful": True,
                        "uniprot_accession" : entry_data["uniprot_accession"],
                        "ena_accession": entry_data["ena_accession"],
                        "uniprot_aa_sequence": entry_data["uniprot_aa_sequence"],
                        }
            else:
                logging.info(
                    f"Error in handling UniProt JSON response for {code}. "
                    + "The structure of the JSON response may have changed. "
                    + "Check handle_uniprot_entry_json function in sd2ba_functions.handle_data"
                    )

                return {"successful": False}

        except:
            logging.info(
                    f"Error in handling UniProt JSON response for {code}. "
                    + "The structure of the JSON response may have changed. "
                    + "Check handle_uniprot_entry_json function in sd2ba_functions.handle_data"
                    )
            return {"successful": False}

    else:
        logging.info(
                f"UniProt request for {code} entry was unsuccessful. "
                + f"Url used was {url}."
                )
        return {"successful": False}


def get_pdb_file(code):
    """
    Get PDB file data for a given code.

    Parameters:
        code (str): The PDB code to retrieve the file data for.

    Returns:
        dict: A dictionary containing the result of the operation. The dictionary has the following keys:
            - "success" (bool): Indicates whether the operation was successful.
            - "url" (str): The URL used to fetch the PDB file data.
            - "data" (requests.Response): The response object containing the PDB file data.

    """

    url = RSCB_PDB_API_URL.format(code=code)

    response = make_request_with_retries(url)

    if response.status_code == 200:
        logging.debug(
                f"RSCB PDB request for {code} pdb file was successful. "
                + f"Url used was {url}."
                )
        return response.text
    else:
        logging.critical(
                f"RSCB PDB request for {code} pdb file was unsuccessful. "
                + f"Url used was {url}. The script has stopped."
                )
        sys.exit("\n** CRITICAL: Request to RSCB PDB failed. Check log file for detailed info. **")

