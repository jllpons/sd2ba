#!/usr/bin/env python3

"""
sd2ba: Segmented Domain to Breakpoint Analysis
"""

import argparse
import logging
import os
import sys


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
            help = 'Specify the logging level. Options: DEBUG, INFO, WARNING, ERROR, CRITICAL. [Default: INFO]'
            )
    parser.add_argument(
            '-v', '--version',
            action='version',
            version=f'sd2ba.py v{__version__}'
            )

    args = parser.parse_args()
    if not args.pdb:
        parser.print_help()
        sys.exit('\n** The pdb, uniprot, codes are required **')

    # Output directory
    output_path = args.o
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = os.path.abspath(output_path) + "/"

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

    pdb_code = args.pdb.upper()

if __name__ == "__main__":
    main()
