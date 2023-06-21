#!/usr/bin/env python3

"""
sd2ba: Segmented Domain to Breakpoint Analysis
"""

import argparse
import logging
import os
import sys

from sd2ba_functions.fetch_data import get_pdb_file


__version__ = '0.0'


def main():

    description = """sd2ba: Segmented Domain to Breakpoint Analysis"""

    parser = argparse.ArgumentParser(
                        description=description,
                        usage="sd2ba.py [options] -pdb STR"
                        )

    parser_codes = parser.add_argument_group('required codes')
    parser_codes.add_argument(
                    '-pdb',
                    metavar='STR',
                    type=str,
                    help='PDB code of the reference protein including the chain',
                    )

    parser_output = parser.add_argument_group('output options')
    parser_output.add_argument(
                    '-o',
                    metavar='STR',
                    type=str,
                    default='sd2ba_output',
                    help = 'path/to/output_directory [Default: s2ba_output]'
                    )

    parser.add_argument(
                    '-l', '--logging',
                    metavar='STR',
                    type=str,
                    default='INFO',
                    help = 'Specify the logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL. [Default: INFO]'
                    )
    parser.add_argument(
                    '-v', '--version',
                    action='version',
                    version=f'sd2ba.py v{__version__}'
            )

    args = parser.parse_args()

    # Check if all required args are present
    if not args.pdb:
        parser.print_help()
        sys.exit('\n** The pdb, uniprot, codes are required **')

    # Output directory
    output_path = args.o
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = os.path.abspath(output_path)

    # Set up logger
    log_level = getattr(logging, args.logging.upper())
    logging.basicConfig(
                filename=f'{output_path}/sd2ba.log',
                encoding='utf-8',
                format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                level=log_level,
                )

    logging.info(f'sd2ba.py v{__version__}')
    logging.info('CMD: ' + ' '.join(sys.argv))
    logging.debug(f'Output directory was created: {output_path}')

    pdb_code = args.pdb.upper()

    # Requesting PDB file to RSCB-PDB and handleling a possible unsuccessful request.
    pdb_response = get_pdb_file(pdb_code)
    if pdb_response['success'] is True:
        logging.debug(
                f'RSCB PDB request for {pdb_code} pdb file was unsuccessful. '
                + f'Url used was {pdb_response["url"]}.'
                )
    else:
        logging.critical(
                f'RSCB PDB request for {pdb_code} pdb file was unsuccessful. '
                + f'Url used was {pdb_response["url"]}. The script has stopped.'
                )
        sys.exit('\n** CRITICAL: Request to RSCB PDB failed. Check log file for detailed info. **')

    # Saving the PDB file to the output directory
    pdb_file_path = output_path + f'/{pdb_code}.pdb'
    with open(pdb_file_path, 'w') as handle:
        handle.write(pdb_response['data'].text)
    logging.debug(f'PDB file for {pdb_code} was saved in {pdb_file_path}')

if __name__ == "__main__":
    main()
