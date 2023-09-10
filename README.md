# sd2ba: Segmented Domain to Breakpoint Analysis

NOTE: this script its currently under development.

---

sd2ba.py is a script that tries to simplify the study of proteins containing various domains that are not contiguous in their amino acid sequence. 

sd2ba.py takes:

1. The **PDB code** of a solved structure from a protein containing various domains that are not contiguous in the sequence.
2. **Two fasta files**, one containing the nucleotide and the other containing
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

### 3rd party programs

- [HMMER](<https://anaconda.org/bioconda/hmmer>): perform a MSA of the input
  proteins using a HMM

```shell
conda install -c bioconda hmmer
```

- [pal2nal](<https://anaconda.org/bioconda/pal2nal>): convert a MSA and the corresponding DNA sequences into a codon alignmet

```shell
conda install -c bioconda pal2nal
```

- [HYPHY](<https://anaconda.org/bioconda/hyphy/>): to analyze a codon aligment with GARD (Genetic Algorithm for Recombination Detection)

```shell
conda install -c bioconda hyphy
```

