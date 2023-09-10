# sd2ba: Segmented Domain to Breakpoint Analysis

NOTE: this script its currently under development.

---

sd2ba.py is a script that tries to simplify the study of proteins containing various domains that are not contiguous in their amino acid sequence. 

sd2ba.py takes:

1. The **PDB code** of a solved structure from a protein containing various domains that are not contiguous in the sequence.
2. The **PFAM code** for the PFAM protein family of the input PDB code protein.
3. **Two fasta files**, one containing the nucleotide and the other containing
   the amino acid sequences, of *n* proteins that share the same domain disposition.

Given this input, sd2ba.py generates a report where the positions of the segmented domains are compared to the positions of  the predicted recombination breakpoints by [GARD](<https://doi.org/10.1093/bioinformatics/btl474>) using a codon aligment from the two fasta files.

## Dependencies

Python (>= 3.9)

The most convenient way of running sd2ba.py is through a conda environment:

```shell
conda create --name sd2ba python=3.9
```

### 3rd party Python libraries:

- [requests](<https://pypi.org/project/requests/>): fetch data from databases
- [biopython](<https://pypi.org/project/biopython/>): handle PDB files

```shell
conda install -c anaconda requests
conda install -c conda-forge biopython
```

### 3rd party programs:

- [HMMER](<https://anaconda.org/bioconda/hmmer>): perform a MSA of the input
  proteins using a HMM
- [pal2nal](<https://anaconda.org/bioconda/pal2nal>): convert a MSA and the corresponding DNA sequences into a codon alignmet
- [HYPHY](<https://anaconda.org/bioconda/hyphy/>): analyze a codon aligment with GARD (Genetic Algorithm for Recombination Detection)

```shell
conda install -c bioconda hmmer
conda install -c bioconda pal2nal
conda install -c bioconda hyphy
```

## Usage

```text
./sd2ba/sd2ba.py -h
usage: sd2ba.py PDB PFAM amino.fasta nucleotide.fasta [options]

sd2ba.py is a tool to identify the breakpoint of a segmented domain in a protein family

positional arguments:
  PDB                   PDB code of the reference protein including the chain
  PFAM                  PFAM code of the domain of interest. The associated HMM will be downloaded.
  amino.fasta           Fasta file with the amino acid sequences of the proteins of interest.
  nucleotide.fasta      Fasta file with the nucleotide sequences of the proteins of interest.

optional arguments:
  -h, --help            show this help message and exit
  -l STR, --logging STR
                        Set the logging level [Default: INFO] [Choices: DEBUG, INFO, WARNING, ERROR, CRITICAL]
  -V, --version         show program's version number and exit

output options:
  -o STR, --output-directory STR
                        Output directory [Default: $CWD/s2ba_output]
```

### Example:

```shell
python sd2ba/sd2ba.py 2QGUA PF05494 test/data/3_aa_sequences.fasta test/data/3_nucleotide_sequences.fasta -o 3_protein_test
```

## Documentation

WIP.
