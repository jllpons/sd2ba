#!/usr/bin/env python3

HMMER_ARGS = r"hmmalign --outformat afa --trim {hmm_filepath} {aa_seq_filepath}"

SED_ARGS_AFA_TO_FASTA = ["sed", r"s/\./-/g"]
