import numpy as np

def LocalAlignment(seq1, seq2, score):
    """
    Calculates the local alignment of two sequences using the Smith-Waterman algorithm.

    Parameters:
    seq1 (str): The first sequence.
    seq2 (str): The second sequence.

    Returns:
    tuple: Aligned sequences (sequence1, sequence2).
    """

    match, mismatch, gap = score  # Unpack the scoring system

    len1, len2 = len(seq1), len(seq2)

    # Initialize score matrix S
    S = np.zeros((len1 + 1, len2 + 1))

    # Fill the score matrix
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if seq1[i - 1] == seq2[j - 1]:
                score = match
            else:
                score = mismatch
            S[i, j] = max(0,
                          (S[i-1, j] - abs(gap)),
                          (S[i, j-1] - abs(gap)),
                          (S[i-1, j-1] + score))

    print("Score matrix:")
    print(S)

    sequence1 = ''
    sequence2 = ''
    [r, s] = np.unravel_index(np.argmax(S), S.shape)

    # Trace back to find the best local alignment
    while r != 0 and s != 0:
        if S[r, s] == 0:
            break
        options = [0,
                  S[r - 1, s] - abs(gap),
                  S[r, s - 1] - abs(gap),
                  S[r - 1, s - 1] + (match if seq1[r-1] == seq2[s-1] else mismatch)]
        index = options.index(max(options))
        if index == 1:
            sequence1 = seq1[r - 1] + sequence1
            sequence2 = "-" + sequence2
            r -= 1
        elif index == 2:
            sequence1 = "-" + sequence1
            sequence2 = seq2[s - 1] + sequence2
            s -= 1
        else:
            sequence1 = seq1[r - 1] + sequence1
            sequence2 = seq2[s - 1] + sequence2
            r -= 1
            s -= 1

    return sequence1, sequence2


# Test function call
print(LocalAlignment(seq2='AATGC', seq1='AGC', score=[2,-1,-2]))

