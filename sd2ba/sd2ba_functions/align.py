#!/usr/bin/env python3

def start_pos_in_trimmed_sequence(complete, trimmed, n=6):
    """
    Returns the position of the start of a trimmed sequence compared to the
    complete sequence.

    Args:
        complete (str): The complete sequence.
        trimmed (str): The trimmed sequence.
        n (int): The number of aminoacids to compare. Default is 6.

    Returns:
        int: The position of the start of the trimmed sequence compared to the
             complete sequence.

    """

    start = n - n
    end = n

    chunk = trimmed[start:end]

    # Don't want gaps in our chunk
    while "-" in chunk:
        start += 1
        end += 1
        chunk = trimmed[start:end]

    trimmed_start = start

    start = n - n
    end = n

    while chunk != complete[start:end]:
        start += 1
        end += 1

    complete_match = start

    return complete_match - trimmed_start


# When did I write this?
# def insert_insertions(protein_with_insertions, protein_without_insertions):
# 
#     n_insertions = protein_with_insertions.count("-")
# 
#     while n_insertions > 0:
#         for index, value in enumerate(protein_with_insertions):
#             if value == "-":
#                 protein_without_insertions = protein_without_insertions[:index] + "-" + protein_without_insertions[index:]
#                 n_insertions -= 1
#                 break
# 
#     return protein_without_insertions
