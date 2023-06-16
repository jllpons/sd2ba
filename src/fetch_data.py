#!/usr/bin/env python3

from api_urls import UNIPROT_API_URL

def make_uniprot_api_call(code):

    url = UNIPROT_API_URL.format(code=code)
    # Make request to the dynamic URL
    # ...

# make_api_call("12345")  # Replace "12345" with the actual dynamic value

