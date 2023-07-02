#!/usr/bin/env python3

import logging
import sys

import requests

from requests.adapters import HTTPAdapter, Retry

from sd2ba_functions.api_urls import RSCB_PDB_API_URL


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
                        allowed_methods=['GET', 'POST'],
                        backoff_factor=delay,
                        )

    adapter = HTTPAdapter(max_retries=retry_strategy)
    http = requests.Session()
    http.mount('https://', adapter)

    response = http.get(url)
    return response


def get_pdb_file(code):
    """
    Get PDB file data for a given code.

    Parameters:
        code (str): The PDB code to retrieve the file data for.

    Returns:
        dict: A dictionary containing the result of the operation. The dictionary has the following keys:
            - 'success' (bool): Indicates whether the operation was successful.
            - 'url' (str): The URL used to fetch the PDB file data.
            - 'data' (requests.Response): The response object containing the PDB file data.

    """

    url = RSCB_PDB_API_URL.format(code=code)

    response = make_request_with_retries(url)

    if response.status_code == 200:
        logging.debug(
                f'RSCB PDB request for {code} pdb file was successful. '
                + f'Url used was {url}.'
                )
        return response
    else:
        logging.critical(
                f'RSCB PDB request for {code} pdb file was unsuccessful. '
                + f'Url used was {url}. The script has stopped.'
                )
        sys.exit('\n** CRITICAL: Request to RSCB PDB failed. Check log file for detailed info. **')

