import numpy as np

def GlobalAlignment(seq1, seq2, score):
    """
    Calculates the global alignment of two sequences using the Needleman-Wunsch algorithm.

    Parameters:
    seq1 (str): The first sequence.
    seq2 (str): The second sequence.
    score (tuple): Scoring system in the form (match, mismatch, gap).

    Returns:
    tuple: Score matrix and aligned sequences (sequence1, sequence2).
    """

    match, mismatch, gap = score  # Unpack the scoring system

    len1, len2 = len(seq1), len(seq2)

    # Initialize score matrix S with penalties for gaps
    S = np.zeros((len1 + 1, len2 + 1))
    S[0, :] = np.arange(0, (len2 + 1) * gap, gap)
    S[:, 0] = np.arange(0, (len1 + 1) * gap, gap)

    # Fill the score matrix
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            a = match if seq1[i - 1] == seq2[j - 1] else mismatch
            S[i, j] = max(
                S[i - 1, j] - abs(gap),  # Gap in seq2
                S[i, j - 1] - abs(gap),  # Gap in seq1
                S[i - 1, j - 1] + a      # Alignment of characters
            )

    print("Score matrix:")
    print(S)

    # Backtracking to find the aligned sequences
    sequence1, sequence2 = '', ''
    r, s = len1, len2

    while r > 0 or s > 0:
        if r == 0:  # Only seq2 continues
            sequence1 = "-" + sequence1
            sequence2 = seq2[s - 1] + sequence2
            s -= 1
        elif s == 0:  # Only seq1 continues
            sequence1 = seq1[r - 1] + sequence1
            sequence2 = "-" + sequence2
            r -= 1
        else:
            # Choose the best direction to go back
            options = [
                (S[r - 1, s] - abs(gap), r - 1, s),  # Gap in seq2
                (S[r, s - 1] - abs(gap), r, s - 1),  # Gap in seq1
                (S[r - 1, s - 1] + (match if seq1[r - 1] == seq2[s - 1] else mismatch), r - 1, s - 1)  # Alignment
            ]
            max_index = max(options, key=lambda x: x[0])

            if max_index[1] == r - 1 and max_index[2] == s:
                sequence1 = seq1[r - 1] + sequence1
                sequence2 = "-" + sequence2
                r -= 1
            elif max_index[1] == r and max_index[2] == s - 1:
                sequence1 = "-" + sequence1
                sequence2 = seq2[s - 1] + sequence2
                s -= 1
            else:
                sequence1 = seq1[r - 1] + sequence1
                sequence2 = seq2[s - 1] + sequence2
                r -= 1
                s -= 1

    print("Aligned sequences:")
    print(sequence1, sequence2)

    return S, sequence1, sequence2


# Test function call
GlobalAlignment('ACGTA', 'TGGGACC', (2, -1, -2))


