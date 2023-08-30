#!/usr/bin/env python3

UNIPROT_UNIREF_JSON_API_URL = "https://rest.uniprot.org/uniref/{code}.json"

UNIPROT_ENTRY_JSON_API_URL = "https://rest.uniprot.org/uniprotkb/{code}.json"

RSCB_PDB_API_URL = "https://files.rcsb.org/view/{code}.pdb"

ENA_NUCLEOTIDE_SEQUENCE_API_URL = "https://www.ebi.ac.uk/ena/browser/api/fasta/{code}"

PFAM_HMM_API_URL = "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{code}?annotation=hmm"

RSCB_PDB_FASTA_API_URL = "https://www.rcsb.org/fasta/entry/{code}/display"

PDB_TO_UNIPROT_JSON_API_URL = "https://rest.uniprot.org/uniprotkb/search?format=json&query=%28{code}%29&size=500"

EBI_MAPPING_PDBE_API_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/{code}"
