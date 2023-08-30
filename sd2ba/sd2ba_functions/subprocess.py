#!/usr/bin/env python3

import logging
import subprocess
import sys


logger = logging.getLogger(__name__)


def run_subprocess(process_name, args, input=None):
    """Run a subprocess and return stdout as a string.

    Args:
        process_name (str): Name of the process to be run.
        args (list): List of arguments to be passed to the process.
        input (str): Input to be passed to the process. Defaults is None.

    Returns:
        stdout (bytes): stdout of the process.

    Raises:
        SystemExit: If the process returns a non-zero exit code.
    """

    if input:
        call = subprocess.run(
                            args,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            input=input,
                            )
    else:
        call = subprocess.run(
                            args,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            )

    logging.info(f"{process_name} call: {call.args}")
    logging.debug(f"{process_name} stdout: {call.stdout.decode('utf-8')}")
    logging.debug(f"{process_name} stderr: {call.stderr.decode('utf-8')}")
    logging.debug(f"{process_name} returncode: {call.returncode}")

    if call.returncode != 0:
        logging.critical(f"{process_name} failed with returncode {call.returncode}")
        logging.critical("script ended")
        sys.exit(f"\n** {process_name} call failed, check log file for more info **")

    return call.stdout

