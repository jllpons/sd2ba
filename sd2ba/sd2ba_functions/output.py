#!/usr/bin/env python3

"""
"""

def write_fasta(path, content):
    """
    """
    with open(path, 'w') as handle:
        handle.write(content)

def generate_representation(sequence, fragments, header=None):
    """
    Generates a representation of the sequence with the fragments. Can be either
    a segmented domain or predicted breakpoints. The representation is a string
    where the segments are represented by two colons (::) and the rest of the
    sequence is represented by the sequence itself.

    Args:
        sequence (str): The sequence to be represented.
        fragments (dict): A dictionary with the fragments. The keys are the
                          fragment ids and the values are dictionaries with the 
                          keys "start" and "end" that represent the start and 
                          end positions of the fragment.
        header (str): The header of the sequence. Default is None.

    Returns:
        str: The representation of the sequence with the fragments.
    """

    segment_positions = []
    for i in fragments:
        # NOTE: Why I'm not adding the end position?
        # lets suppose we have two segments
        #     segment 1: start = 1, end = 5
        #     segment 2: start = 6, end = 10
        # segment 1 containts 1, 2, 3, 4, 5
        # segment 2 containts 6, 7, 8, 9, 10
        # We don't want a representation like this
        #     ::123456::::78910::
        # We want a representation like this
        #     ::12345::678910::

        if "aminoacid_rounded" in list(fragments[i].keys()):
            start = fragments[i]["aminoacid_rounded"]["start"]
            segment_positions.append(start)

        else:
            start = fragments[i]["start"]
            segment_positions.append(start)

    for i in fragments:
        # NOTE: But at the same time:
        # lets suppose we have two segments
        #     segment 1: start = 1, end = 5
        #     segment 2: start = 7, end = 10
        # segment 1 containts 1, 2, 3, 4, 5
        # segment 2 containts 7, 8, 9, 10
        # We don't want a representation like this
        #     ::123456::78910::
        # We want a representation like this
        #     ::12345::6::78910::

        if "aminoacid_rounded" in fragments[i].keys():
            end = fragments[i]["aminoacid_rounded"]["end"]
            if end not in segment_positions and end + 1 not in segment_positions:
                segment_positions.append(end+1)

        else:
            end = fragments[i]["end"]
            if end + 1 not in segment_positions:
                segment_positions.append(end-1)

    representation = []
    accumulator = -1
    for i in range(len(sequence)):
        accumulator += 1

        if accumulator+1 in segment_positions:
            representation.append('::')
            representation.append(sequence[i])
        else:
            representation.append(sequence[i])

    if header:

        return (header
                + "\n"
                + "".join(representation))

    else:

        return "".join(representation)


def generate_text_report(dict_report, pdb_code):
    """
    Generates a text report with the results of the GARD analysis.

    Args:
        dict_report (dict): A dictionary with the results of the GARD analysis.
        pdb_code (str): The PDB code of the solved structure.

    Returns:
        str: The text report.
    """

    first_header = f"1. Domain data for the solved structure: {pdb_code}\n\n"
    solved_residues = (f"\t- {dict_report['domain_data']['n_solved_residues']} "
                       + f"solved residues: [{dict_report['domain_data']['list_solved_residues'][0]}..{dict_report['domain_data']['list_solved_residues'][-1]}]\n")
    domains_header = f"\t- Domains ({len(dict_report['domain_data']['domains'])}):\n"
    segments = [f"\t\t- {dict_report['domain_data']['domains'][i]['id']}: amino [{dict_report['domain_data']['domains'][i]['start']}..{dict_report['domain_data']['domains'][i]['end']}]\n" for i in dict_report["domain_data"]["domains"]]
    domains = domains_header + "".join(segments)
    #first_subheader = f"\n1.1. Domain data representation for {pdb_code}\n\n\t"
    first_representation = f"{dict_report['domain_data_representation']}\n"
    second_header = "\n2. GARD results for using the codon aligment of the input proteins\n\n"
    predicted_breakpoints_header = f"\t- Predicted breakpoints ({len(dict_report['gard_results'])}):\n"
    breakpoints = [f"\t\t- {i}: amino [{dict_report['gard_results'][i]['aminoacid_rounded']['start']}..{dict_report['gard_results'][i]['aminoacid_rounded']['end']}], nucleotide [{dict_report['gard_results'][i]['nucleotide']['start']}..{dict_report['gard_results'][i]['nucleotide']['end']}]\n" for i in dict_report["gard_results"]]
    predicted_breakpoints = predicted_breakpoints_header + "".join(breakpoints)
    second_subheader = f"\n2.1. GARD predicted breakpoints positions when the trimed amino acids in the {pdb_code} sequence are taken into account:\n\n"
    little_explanation = f"\t(This means that the number of residues trimmed by HMMER ({dict_report['start_pos_of_aligned_aa_sequence']}) are added to the predicted breakpoints positions)\n\n"
    breakpoints_w_trimmed = [f"\t\t- {i}: amino [{dict_report['gard_results_with_aligned_start_pos'][i]['aminoacid_rounded']['start']}..{dict_report['gard_results_with_aligned_start_pos'][i]['aminoacid_rounded']['end']}], nucleotide [{dict_report['gard_results_with_aligned_start_pos'][i]['nucleotide']['start']}..{dict_report['gard_results_with_aligned_start_pos'][i]['nucleotide']['end']}]\n" for i in dict_report["gard_results_with_aligned_start_pos"]]
    predicted_breakpoints_w_trimmed = predicted_breakpoints_header + "".join(breakpoints_w_trimmed)
    third_subheader = "\n2.2 Comparison of the representations of the segmented domain and the predicted breakpoints\n\n"
    second_representation = (dict_report["start_pos_of_aligned_aa_sequence"] * "A"
                             + f"{dict_report['gard_results_representation'][pdb_code.upper()]}")
    second_representation_fasta = f">{pdb_code} | UNIPROT SEQUENCE WITH TRIMMED AA | GARD results representation\n{second_representation}\n"
    fourth_subheader = "\n2.3. Representation of the predicted breakpoints for all of the input proteins\n\n"
    representations = [f">{i}\n{dict_report['gard_results_representation'][i]}\n" for i in dict_report["gard_results_representation"]]

    return (first_header
            + solved_residues
            + domains
            #+ first_subheader
            #+ first_representation
            + second_header
            + predicted_breakpoints
            + second_subheader
            + little_explanation
            + predicted_breakpoints_w_trimmed
            + third_subheader
            + first_representation
            + second_representation_fasta
            + fourth_subheader
            + "".join(representations)
            )




