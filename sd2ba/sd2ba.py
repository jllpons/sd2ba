#!/usr/bin/env python3


"""
sd2ba: Segmented Domain to Breakpoint Analysis
"""


import argparse
import json
import logging
import os
import sys

from sd2ba_functions.fetch_data import (
        get_ena_nucleotide_sequence, get_pdb_file, get_hmm_file,
        get_aa_sequence_from_pdb, get_up_code_from_pdb_code, get_uniprot_entry_data,
        get_domain_info_from_pdb,
        )
from sd2ba_functions.handle_data import (
        get_solved_residues_from_pdb, read_single_fasta, read_multiple_fasta,
        handle_gard_json_output, map_domain_from_cath, handle_domain_mapping,
        )
from sd2ba_functions.SCRIPT_ARGS import (
        HMMER_ARGS, SED_ARGS_AFA_TO_FASTA, PAL2NAL_ARGS, HYPHY_GARD_ARGS,
        )
from sd2ba_functions.subprocess import run_subprocess
from sd2ba_functions.identity import calculate_identity_score
from sd2ba_functions.align import start_pos_in_trimmed_sequence
from sd2ba_functions.output import generate_representation, generate_text_report


__version__ = "0.0.0"


class Protein:

    def __init__(self, header: str):
        self.header = header
        self.aa_sequence = ""
        self.nt_sequence = ""
        self.aligned_aa_sequence = ""
        self.identity_score = 0


class ReferenceProtein(Protein):

    def __init__(self, header: str):
        super().__init__(header)
        self.aa_sequence_pdb = ""
        self.aa_with_known_positions = {}
        self.aligned_aa_sequence_pdb = ""
        self.domain_maping = {}
        self.start_pos_of_aligned_aa_sequence = 0


def main():

    description = """sd2ba.py is a tool to identify the breakpoint of a segmented domain in a protein family"""

    parser = argparse.ArgumentParser(
                        description=description,
                        usage="sd2ba.py PDB PFAM amino.fasta nucleotide.fasta [options]"
                        )

    parser.add_argument(
            "pdb_code",
            metavar="PDB",
            type=str,
            help="PDB code of the reference protein including the chain",
            )
    parser.add_argument(
            "pfam_code",
            metavar="PFAM",
            type=str,
            help="PFAM code of the domain of interest. The associated HMM will be downloaded.",
            )
    parser.add_argument(
            "aa_fasta",
            metavar="amino.fasta",
            type=str,
            help="Fasta file with the amino acid sequences of the proteins of interest.",
            )
    parser.add_argument(
            "nt_fasta",
            metavar="nucleotide.fasta",
            type=str,
            help="Fasta file with the nucleotide sequences of the proteins of interest.",
            )

    parser_output = parser.add_argument_group("output options")
    parser_output.add_argument(
                    "-o", "--output-directory",
                    metavar="STR",
                    type=str,
                    default="sd2ba_output",
                    help = "Output directory [Default: $CWD/s2ba_output]"
                    )

    parser.add_argument(
                    "-l", "--logging",
                    metavar="STR",
                    type=str,
                    default="INFO",
                    help = "Set the logging level [Default: INFO] [Choices: DEBUG, INFO, WARNING, ERROR, CRITICAL]"
                    )
    parser.add_argument(
                    "-V", "--version",
                    action="version",
                    version=f"%(prog)s v{__version__}",
            )

    args = parser.parse_args()

    # PDB code
    pdb_code = args.pdb_code.upper()
    if len(pdb_code) != 5:
        sys.exit(f"\n** The PDB code of the reference protein must specify the chain and provided was: {pdb_code} **")
    # PFAM code
    pfam_code = args.pfam_code.upper()
    # Amino fasta file
    aa_fasta = args.aa_fasta
    if not os.path.isfile(aa_fasta):
        sys.exit(f"\n** The provided amino fasta file does not exist: {aa_fasta} **")
    # Nucleotide fasta file
    nt_fasta = args.nt_fasta
    if not os.path.isfile(nt_fasta):
        sys.exit(f"\n** The provided nucleotide fasta file does not exist: {nt_fasta} **")

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

    logging.info(f"{__file__} v{__version__}")
    logging.info("script started")
    logging.info("CMD: " + " ".join(sys.argv))
    logging.debug(f"Output directory was created: {output_path}")

    # Requesting PDB file to RSCB-PDB.
    # Any unsuccessful request will result in critical error and the exit of the script.
    # The PDB file is saved in the output directory
    pdb_file = output_path + f"/{pdb_code}.pdb"
    with open(pdb_file, "w") as handle:
        # get_pdb_file() returns a string with the obtained PDB file
        # If the request was unsuccessful, the script will exit
        handle.write(get_pdb_file(pdb_code[:4]))
    logging.debug(f"PDB file for {pdb_code} was saved in {pdb_file}")

    ref_protein = ReferenceProtein(header=pdb_code)
    try:
        # get_solved_residues_from_pdb() accesses the PDB file
        # and returns a dictionary with the solved residues
        ref_protein.aa_with_known_positions=get_solved_residues_from_pdb(pdb_code, pdb_file)
    except:
        logging.critical(f"Error while parsing the PDB file of {pdb_code}")
        logging.critical("script ended")
        sys.exit("\n** Error while parsing the PDB file, check the log file for more information **")

    ref_protein.aa_sequence_pdb = read_single_fasta(
            get_aa_sequence_from_pdb(pdb_code[:4])
            )["sequence"]
    uniprot_entry = get_uniprot_entry_data(get_up_code_from_pdb_code(pdb_code[:4]))
    if uniprot_entry["successful"] == False:
        logging.critical(f"Error accessing the Uniprot entry of {pdb_code}")
        logging.critical("script ended")
        sys.exit("\n** Error accessing the Uniprot entry, check the log file for more information **")
    ref_protein.aa_sequence = uniprot_entry["uniprot_aa_sequence"]

    ena_entry = get_ena_nucleotide_sequence(uniprot_entry["ena_accession"])
    if ena_entry["successful"] == False:
        logging.critical(f"Error accessing the ENA entry of {pdb_code}")
        logging.critical("script ended")
        sys.exit("\n** Error accessing the ENA entry, check the log file for more information **")
    ref_protein.nt_sequence = ena_entry["fasta"]["sequence"]

    logging.info(f"Reference protein {pdb_code} has {len(ref_protein.aa_with_known_positions)} solved residues: "
            + f"[{list(ref_protein.aa_with_known_positions.keys())[0]}"
            + f"..{list(ref_protein.aa_with_known_positions.keys())[-1]}]"
            )

    hmm_file = output_path + f"/{pfam_code}.hmm"
    with open(hmm_file, "w") as handle:
        # get_hmm_file() returns a string with the obtained HMM file
        # If the request was unsuccessful, the script will exit
        handle.write(get_hmm_file(pfam_code))
    logging.debug(f"HMM file for {pfam_code} was saved in {hmm_file}")

    # Parsing the amino and nucleotide fasta files
    aa_sequences = {}
    with open(aa_fasta, "r") as handle:
        records = read_multiple_fasta(handle.read())

        if records["successful"] == False:
            logging.critical(f"Error while parsing the amino fasta file {aa_fasta}")
            logging.critical("script ended")
            sys.exit("\n** Error while parsing the amino fasta file, check the log file for more information **")

        aa_sequences = records["data"]

    logging.info(f"{len(aa_sequences)} amino acid sequences were read from {aa_fasta}")

    nt_sequences = {}
    with open(nt_fasta, "r") as handle:
        records = read_multiple_fasta(handle.read())

        if records["successful"] == False:
            logging.critical(f"Error while parsing the nucleotide fasta file {nt_fasta}")
            logging.critical("script ended")
            sys.exit("\n** Error while parsing the nucleotide fasta file, check the log file for more information **")

        nt_sequences = records["data"]

    logging.info(f"{len(nt_sequences)} nucleotide sequences were read from {nt_fasta}")

    # The number of amino acid and nucleotide sequences must be the same
    if len(aa_sequences) != len(nt_sequences):
        logging.critical("The number of amino acid and nucleotide sequences must be the same")
        logging.critical("script ended")
        sys.exit("\n** The number of amino acid and nucleotide sequences must be the same **")

    # Loading the amino acid and nucleotide sequences into Protein objects
    proteins = []
    for header in aa_sequences.keys():
        p = Protein(header=header)
        p.aa_sequence = aa_sequences[header]
        p.nt_sequence = nt_sequences[header]
        proteins.append(p)

    # The reference protein will be the first one in the list
    proteins.insert(0, ref_protein)
    logging.info(f"Reference protein {pdb_code} was inserted as the first protein in the list")

    # At the end, we'll insert the reference protein again, but with the PDB amino acid sequence.
    # It won't be included in the codon aligment nor in the GARD analysis.
    # But we'll use it to understand how the known domain positions are distributed in the MSA
    proteins.append(Protein(header=f"{ref_protein.header}_PDBsequence"))
    proteins[-1].aa_sequence=ref_protein.aa_sequence_pdb

    logging.info(f"Reference protein {pdb_code} PDB sequence was inserted as the last protein in the list")

    input_with_pdb_aa_fasta = output_path + f"/input_with_{pdb_code}_amino.fasta"
    input_with_pdb_nt_fasta = output_path + f"/input_with_{pdb_code}_nucleotide.fasta"

    # Rerwiting the input files with the reference protein as the first one
    with open(input_with_pdb_aa_fasta, "w") as handle:
        for p in proteins:
            handle.write(f">{p.header}\n{p.aa_sequence}\n")

    # The last protein in the list is the reference protein with the PDB amino acid sequence
    # It doesn't have a nucleotide sequence, so we'll remove it from the list
    proteins.pop(-1)

    with open(input_with_pdb_nt_fasta, "w") as handle:
        for p in proteins:
            handle.write(f">{p.header}\n{p.nt_sequence}\n")

    logging.debug(
            f"Nuclotide and amino acid sequences file were rewritten with {pdb_code} as the first sequence "
            + f"and saved in {input_with_pdb_aa_fasta} and {input_with_pdb_nt_fasta}"
            )

    # Running HMMER to make a multiple sequence alignment
    hmmalign_stdout = run_subprocess(
                        "HMMALIGN",
                        HMMER_ARGS
                            .format(
                                hmm_filepath=hmm_file,
                                aa_seq_filepath=input_with_pdb_aa_fasta,
                                )
                            .split(),
                            )
    logging.info("HMMER's HMMALIGN  was used to make a multiple sequence alignment")

    # The hmmalign output is written in afa format. We will convert it to fasta
    # and save it to a file
    hmmalign_output = output_path + "/hmmalign_output.fasta"
    with open(hmmalign_output, "w") as handle:
        sed_stdout = run_subprocess(
                        "SED",
                        SED_ARGS_AFA_TO_FASTA,
                        input=hmmalign_stdout,
                        )
        handle.write(sed_stdout.decode("utf-8"))

    logging.info("SED was called to convert the hmmalign output to fasta format")
    logging.debug(f"HMMER output was saved in {hmmalign_output}")

    # Loading the multiple sequence alignment sequences into the Protein objects
    msa_proteins = []
    with open(hmmalign_output, "r") as handle:
        records = read_multiple_fasta(handle.read())

        if records["successful"] == False:
            logging.critical(f"Error while parsing the multiple sequence alignment file {hmmalign_output}")
            logging.critical("script ended")
            sys.exit("\n** Error while parsing the multiple sequence alignment file, check the log file for more information **")

        msa_proteins = records["data"]

    logging.info(f"{len(msa_proteins)} sequences were read from {hmmalign_output}")

    for p in proteins:
        p.aligned_aa_sequence = msa_proteins[p.header]
        p.identity_score = calculate_identity_score(
                            query=ref_protein.aligned_aa_sequence,
                            subject=p.aligned_aa_sequence,
                            )

    # Identifying the reference protein's PDB aa sequence and saving it in the Protein object
    for p in proteins:
        if p.header == f"{ref_protein.header}_PDBsequence":
            ref_protein.aligned_aa_sequence = p.aa_sequence
            proteins.remove(p)
            logging.info(f"Reference protein {pdb_code} PDB sequence was identified, saved and removed from the aligment.")
            break

    proteins.sort(key=lambda x: x.identity_score, reverse=True)

    logging.debug("Identity scores were calculated for each sequence and proteins were sorted by that value")
    logging.debug(
                "Identity scores: \n"
                + "\n".join(
                    [i.header + ": " + str(i.identity_score) for i in proteins]
                    )
                )
    # NOTE: This is the MSA that will be piped into PAL2NAL
    msa_sorted_by_identity = output_path + "/hmmalign_output_sorted_by_identity.fasta"
    with open(msa_sorted_by_identity, "w") as handle:
        for p in proteins:
            handle.write(f">{p.header}\n{p.aligned_aa_sequence}\n")

    logging.info(f"MSA was sorted by identity score and saved in {msa_sorted_by_identity}")

    # pal2nal will be used to convert the multiple sequence alignment into a codon alignment
    pal2nal_stdout = run_subprocess(
                        "PAL2NAL",
                        PAL2NAL_ARGS
                            .format(
                                msa_aa_filepath=msa_sorted_by_identity,
                                nt_filepath=input_with_pdb_nt_fasta,
                                )
                            .split(),
                            )
    logging.info("PAL2NAL was used to convert the multiple sequence alignment into a codon alignment")

    pal2nal_output = output_path + "/pal2nal_output.fasta"
    with open(pal2nal_output, "w") as handle:
        handle.write(pal2nal_stdout.decode("utf-8"))

    logging.debug(f"PAL2NAL output was saved in {pal2nal_output}")

    # Running hyphy's GARD to detect recombination breakpoints
    print("\n** INFO: Running HYPHY-GARD to detect recombination breakpoints. It may take a while **\n")
    gard_stdout = run_subprocess(
                    "HYPHY-GARD",
                    HYPHY_GARD_ARGS
                        .format(
                            pal2nal_filepath=pal2nal_output,
                            )
                        .split(),
                        )
    logging.info("HYPHY-GARD was used to detect recombination breakpoints")
    logging.debug("HYPHY-GARD stdoutput: \n" + gard_stdout.decode("utf-8"))

    # Parsing the HYPHY-GARD output
    HYPHY_GARD_JSON_OUTPUT = output_path + "/pal2nal_output.fasta.GARD.json"
    logging.debug(f"Parsing HYPHY-GARD JSON output from {HYPHY_GARD_JSON_OUTPUT}")

    breakpoint_info = handle_gard_json_output(HYPHY_GARD_JSON_OUTPUT)
    logging.debug(f"Breakpoint info: \n{breakpoint_info}\n")
    # Sequences used in GARD analysis has been trimmed. We need to find the
    # start position of the aligned sequence in the original sequences
    ref_protein.start_pos_of_aligned_aa_sequence = start_pos_in_trimmed_sequence(
                                                    complete=ref_protein.aa_sequence,
                                                    trimmed=ref_protein.aligned_aa_sequence,
                                                    )
    logging.info(f"Start position of the aligned sequence in the original sequence was identified for {pdb_code}")
    breakpoint_info_with_aligned_start_pos = handle_gard_json_output(
                                                HYPHY_GARD_JSON_OUTPUT,
                                                start_position=ref_protein.start_pos_of_aligned_aa_sequence,
                                                )
    logging.debug(f"Breakpoint info with aligned start pos: \n{breakpoint_info_with_aligned_start_pos}\n")

    cath_domain_data = get_domain_info_from_pdb(pdb_code[:4])
    logging.info(f"Domain info was retrieved from CATH for {pdb_code}")
    logging.debug(f"Domain info: \n{cath_domain_data}\n")

    ref_protein.domain_mapping = map_domain_from_cath(cath_domain_data)

    # Crafting the report
    report = {
            "domain_data" : {
                "n_solved_residues" : len(ref_protein.aa_with_known_positions),
                "list_solved_residues" : [i for i in ref_protein.aa_with_known_positions],
                "domains" : handle_domain_mapping(ref_protein.domain_mapping),
                },
            "domain_data_representation" : generate_representation(
                                                sequence=ref_protein.aa_sequence_pdb,
                                                fragments=handle_domain_mapping(ref_protein.domain_mapping),
                                                start_pos=[i for i in ref_protein.aa_with_known_positions][0],
                                                ),
            "gard_results" : breakpoint_info["data"],
            "start_pos_of_aligned_aa_sequence" : ref_protein.start_pos_of_aligned_aa_sequence,
            "gard_results_with_aligned_start_pos" : breakpoint_info_with_aligned_start_pos["data"],
            "gard_results_representation" : {
                i.header : generate_representation(
                    i.aligned_aa_sequence,
                    breakpoint_info_with_aligned_start_pos["data"])
                for i in proteins
                },
            }

    report_json_path = output_path + "/report.json"
    with open(report_json_path, "w") as handle:

        json.dump(report, handle, indent=4)

    report_txt_path = output_path + "/report.txt"
    with open(report_txt_path, "w") as handle:
        handle.write(generate_text_report(report, pdb_code))

    logging.info(f"Report was generated and saved in {report_json_path} and {report_txt_path}")
    logging.info(f"Program finished successfully")


if __name__ == "__main__":
    main()
