import numpy as np


def BLOSUM_calculation(sequences):
    """
    Function to calculate the BLOSUM matrix based on a given set of sequences.

    Parameters:
    sequences (list of str): A list of sequences (strings) for which to calculate the BLOSUM matrix.

    Returns:
    np.array: The BLOSUM matrix.
    """

    # Step 1: Get the number of sequences, sequence lengths, and the number of possible pairs
    num_sequences = len(sequences)
    seq_length = len(sequences[0])
    num_pairs = num_sequences * (num_sequences - 1) // 2 * seq_length

    # Step 2: Identify the unique amino acids in the sequences and sort them alphabetically
    amino_acids = sorted(set(''.join(sequences)))

    # Step 3: Initialize a frequency matrix F with zeros
    F = np.zeros((len(amino_acids), len(amino_acids)))

    # Step 5: Compute the frequency matrix F
    split_sequences = [list(seq) for seq in sequences]
    for i in range(0, len(amino_acids)):
        for j in range(0, len(amino_acids)):
            if i == j:
                for k in range(seq_length):
                    ni = sum(1 for m in range(3) if split_sequences[m][k] == amino_acids[i])
                    F[i, i] += ni * (ni - 1) / 2
            else:
                for k in range(seq_length):
                    ni = sum(1 for m in range(len(amino_acids)) if split_sequences[m][k] == amino_acids[i])
                    nj = sum(1 for m in range(len(amino_acids)) if split_sequences[m][k] == amino_acids[j])
                    F[i, j] += ni * nj
    print("Frequencies of amino acid pairs are:\n", F)

    # Step 6: Normalize the frequency matrix (divide by total pairs)
    Q = F / num_pairs
    print("Normovaná matice četnosti znaků je:\n", Q)

    # Step 7: Compute marginal sums
    P = np.zeros((len(Q)))
    for i in range(len(Q)):
        P[i] = Q[i, i] + sum(Q[i, j] / 2 for j in range(len(Q)) if i != j)

    print("Marginální součty jsou:\n", P)


    # Step 8: Compute the BLOSUM matrix
    B = np.zeros((Q.shape))
    for i in range(len(amino_acids)):
        for j in range(len(amino_acids)):
            if i == j:
                if Q[i, j] > 0:
                    B[i, j] = 2* np.log2(Q[i,j]/P[i]**2)
                else:
                    B[i, j] = 0
            else:
                if Q[i, j] > 0:
                    B[i, j] = 2*np.log2(Q[i,j]/(2*P[i]*P[j]))
                else:
                    B[i, j] = 0
    print("BLOSUM matice je:\n", np.round(B))
    return B



BLOSUM_calculation(['CAABABA', 'BBABCBB', 'AAABCBA'])