#!/usr/bin/env python3

"""
"""

def write_fasta(path, content):
    """
    """
    with open(path, 'w') as handle:
        handle.write(content)

def generate_representation(sequence, fragments):
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
            if end + 1 not in segment_positions:
                segment_positions.append(end+1)

        else:
            end = fragments[i]["end"]
            if end + 1 not in segment_positions:
                segment_positions.append(end-1)

    representation = []
    accumulator = -1
    for i in range(len(sequence)):
        if sequence[i] != '-':
            accumulator += 1

        if accumulator+1 in segment_positions:
            representation.append('::')
            representation.append(sequence[i])
        else:
            representation.append(sequence[i])

    return "".join(representation)
