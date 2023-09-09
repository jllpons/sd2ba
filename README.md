# sd2ba: Segmented Domain to Breakpoint Analysis

NOTE: this script its currently under development.

---

sd2ba.py is a script that tries to simplify the study of proteins containing various domains that are not contiguous in their amino acid sequence. 

sd2ba.py takes:

1. The **PDB code** of a solved structure from a protein containing various domains that are not contiguous in the sequence.
2. **Two fasta files**, one containing the nucleotide and the other containing
   the amino acid sequences, of *n* proteins that share the same domain disposition.

Given this input, sd2ba.py generates a report where the positions of the segmented domains are compared to the positions of  the predicted recombination breakpoints by [GARD](<https://doi.org/10.1093/bioinformatics/btl474>) using a codon algiment from the two fasta files.

