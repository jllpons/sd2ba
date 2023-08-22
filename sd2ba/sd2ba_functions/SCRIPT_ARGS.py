#!/usr/bin/env python3

HMMER_ARGS = r"hmmalign --outformat afa --trim {hmm_filepath} {aa_seq_filepath}"

SED_ARGS_AFA_TO_FASTA = ["sed", r"s/\./-/g"]

PAL2NAL_ARGS = r"pal2nal.pl {msa_aa_filepath} {nt_filepath} -output fasta -codontable 11"

HYPHY_GARD_ARGS = r"hyphy gard {pal2nal_filepath}"
