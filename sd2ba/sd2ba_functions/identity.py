#!/usr/bin/env python3

# TODO: maybe a matrix could be applied here
def calculate_identity_score(query, subject):
    """
    Calculate the identity score between two sequences.

    Parameters
        query (str): The query sequence.
        subject (str): The subject sequence.

    Returns
        score (int): The identity score.
    """

    # Calculate the identity score
    score = 0
    for i, j in zip(query, subject):
        if i != "-":
            if i == j:
                score += 1

    return score

