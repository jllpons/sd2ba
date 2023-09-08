#!/usr/bin/env python3


import gzip
import logging
import sys

import requests
from requests.adapters import HTTPAdapter, Retry

from sd2ba_functions.API_URL import (
        UNIPROT_UNIREF_JSON_API_URL, RSCB_PDB_API_URL,
        UNIPROT_ENTRY_JSON_API_URL, ENA_NUCLEOTIDE_SEQUENCE_API_URL,
        PFAM_HMM_API_URL, RSCB_PDB_FASTA_API_URL, PDB_TO_UNIPROT_JSON_API_URL,
        EBI_MAPPING_PDBE_API_URL,
        )
from sd2ba_functions.handle_data import handle_uniref_json, handle_uniprot_entry_json, read_single_fasta


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
    Fetches UniRef data in JSON format for a given UniRef code.

    Parameters:
        code (str): The UniRef code for the desired UniRef cluster.

    Returns:
        requests.Response.text: The JSON response text.

    Raises:
        SystemExit: If the request to UniProt for UniRef data fails or if there is an error
                    in handling the UniRef JSON response.

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
                    f"Error in handling UniProt-UniRef JSON response for {code}. "
                    + "The structure of the JSON response could have changed. "
                    + "Check handle_uniref_json() in sd2ba_functions.handle_data"
                    )
            sys.exit("\n** CRITICAL: Error in handling UniProt-UniRef response. Check log file for detailed info. **")

    logging.critical(
            f"UniProt request for {code} was unsuccessful. "
            + f"Url used was {url}. The script has stopped."
            )
    sys.exit("\n** CRITICAL: Request for UniRef data to UniProt failed. Check log file for detailed info. **")


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
        str: The PDB file data in text format.

    Raises:
        SystemExit: If the request to RSCB PDB for the PDB file fails.
    """

    url = RSCB_PDB_API_URL.format(code=code)

    response = make_request_with_retries(url)

    if response.status_code == 200:
        logging.debug(
                f"RSCB PDB request for {code} pdb file was successful. "
                + f"Url used was {url}."
                )
        return response.text

    logging.critical(
            f"RSCB PDB request for {code} pdb file was unsuccessful. "
            + f"Url used was {url}. The script has stopped."
            )
    sys.exit("\n** CRITICAL: Request to RSCB PDB failed. Check log file for detailed info. **")


def get_ena_nucleotide_sequence(code):
    """
    Get ENA nucleotide sequence for a given code.

    Parameters:
        code (str): The ENA code to retrieve the sequence for.

    Returns:
        dict: A dictionary containing the result of the operation. The dictionary has the following keys:
            - "success" (bool): Indicates whether the operation was successful.
            - "fasta" (str): The fasta sequence for the ENA entry.
    """

    url = ENA_NUCLEOTIDE_SEQUENCE_API_URL.format(code=code)

    response = make_request_with_retries(url)

    if response.status_code == 200:
        logging.debug(
                f"ENA request for {code} nucleotide sequence was successful. "
                + f"Url used was {url}."
                )

        ena_fasta = read_single_fasta(response.text)

        if ena_fasta["successful"] == False:
            logging.info(
                    f"FastA sequence for {code} nucleotide sequence was not found. "
                    + f"Url used was {url}."
                    )
            return {"successful": False}

        return {
                "successful": True,
                "fasta": ena_fasta,
                }

    logging.info(
            f"ENA request for {code} nucleotide sequence was unsuccessful. "
            + f"Url used was {url}."
            )
    return {"successful": False}


def get_hmm_file(pfam_code):
    """
    Get HMM file data for a given PFAM code.

    Parameters:
        pfam_code (str): The PFAM code to retrieve the file data for.

    Returns:
        str: The HMM file data in text format.

    Raises:
        SystemExit: If the request to PFAM for the HMM file fails.
    """

    url = PFAM_HMM_API_URL.format(code=pfam_code)

    response = make_request_with_retries(url)

    if response.status_code == 200:
        logging.debug(
                f"PFAM request for {pfam_code} hmm file was successful. "
                + f"Url used was {url}."
                )
        try:
            # Trust issues
            response = gzip.decompress(response.content)
            return response.decode("utf-8")

        except:
            return response.text

    logging.critical(
            f"PFAM request for {pfam_code} hmm file was unsuccessful. "
            + f"Url used was {url}. The script has stopped."
            )
    sys.exit("\n** CRITICAL: Request to PFAM failed. Check log file for detailed info. **")


def get_aa_sequence_from_pdb(pdbcode):
    """
    Get amino acid sequence from RSCB PDB for a given code.

    Parameters:
        pdbcode (str): The PDB code to retrieve the sequence for.

    Returns:
        str: The amino acid sequence in text format.

    Raises:
        SystemExit: If the request to RSCB PDB for the PDB file fails.
    """

    url = RSCB_PDB_FASTA_API_URL.format(code=pdbcode)

    response = make_request_with_retries(url)

    if response.status_code == 200:
        logging.debug(
                f"RSCB PDB request for {pdbcode} fasta amino sequence was successful. "
                + f"Url used was {url}."
                )

        return response.text

    logging.critical(
            f"RSCB PDB request for {pdbcode} fasta amino sequence was unsuccessful. "
            + f"Url used was {url}. The script has stopped."
            )

    sys.exit("\n** CRITICAL: Request to RSCB PDB failed. Check log file for detailed info. **")


def get_up_code_from_pdb_code(pdbcode):
    """
    Get UniProt code from RSCB PDB for a given code.

    Parameters:
        pdbcode (str): The PDB code to retrieve the UniProt code for.

    Returns:
        str: The UniProt code in text format.

    Raises:
        SystemExit: If the request to UniProt for an UniProt code from PDB code fails.
    """

    url = PDB_TO_UNIPROT_JSON_API_URL.format(code=pdbcode)

    response = make_request_with_retries(url)

    if response.status_code == 200:
        logging.debug(
                f"UniProt request for {pdbcode} UniProt code was successful. "
                + f"Url used was {url}."
                )

        results = response.json()

        if not len(results):
            logging.critical(
                    f"UniProt request for {pdbcode} UniProt code was unsuccessful. "
                    + f"Url used was {url}. The script has stopped."
                    )

            sys.exit("\n** CRITICAL: Request to UniProt for a UniProt code "
                    + " from PDB code failed. Check log file for detailed info. **")

        return results["results"][0]["primaryAccession"]

    logging.critical(
            f"RSCB PDB request for {pdbcode} UniProt code was unsuccessful. "
            + f"Url used was {url}. The script has stopped."
            )

    sys.exit("\n** CRITICAL: Request to RSCB PDB failed. Check log file for detailed info. **")


def get_domain_info_from_pdb(pdb_code):
    """
    Get domain info from EBI for a given code.

    Parameters:
        pdbcode (str): The PDB code to retrieve the domain info for.

    Returns:
        dict: The domain info in json format.

    Raises:
        SystemExit: If the request to EBI for domain info fails.
    """

    url = EBI_MAPPING_PDBE_API_URL.format(code=pdb_code)

    response = make_request_with_retries(url)

    if response.status_code == 200:
        logging.debug(
                f"EBI request for {pdb_code} domain info was successful. "
                + f"Url used was {url}."
                )

        results = response.json()

        if not len(results):
            logging.critical(
                    f"EBI request for {pdb_code} domain info was unsuccessful. "
                    + f"Url used was {url}. The script has stopped."
                    )

            sys.exit("\n** CRITICAL: Request to EBI for domain info failed. "
                    + "Check log file for detailed info. **")

        return results[pdb_code.lower()]

    logging.critical(
            f"EBI request for {pdb_code} domain info was unsuccessful. "
            + f"Url used was {url}. The script has stopped."
            )

    sys.exit("\n** CRITICAL: Request to EBI failed. Check log file for detailed info. **")

