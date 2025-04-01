# Bioinformatics Algorithms Repository

This repository contains implementations of fundamental bioinformatics algorithms for sequence alignment and scoring matrix calculation. The included scripts are useful for comparing biological sequences and analyzing sequence similarity.

## 📂 Repository Structure
```
📁 bioinformatics-algorithms/
│-- GlobalAlignment.py        # Needleman-Wunsch algorithm for global sequence alignment
│-- LocalAlignment.py         # Smith-Waterman algorithm for local sequence alignment
│-- BLOSUM_calculation.py     # BLOSUM matrix calculation
│-- README.md                 # Documentation
```

## 📜 Algorithms Implemented

### 1️⃣ Global Alignment (Needleman-Wunsch)
**File:** `GlobalAlignment.py`

This script performs **global sequence alignment** using the Needleman-Wunsch algorithm, which is widely used for aligning two full-length biological sequences.

#### 🔹 Function: `GlobalAlignment(seq1, seq2, score)`
- **Parameters:**
  - `seq1` (str): The first sequence.
  - `seq2` (str): The second sequence.
  - `score` (tuple): Scoring system in the form `(match, mismatch, gap)`.
- **Returns:**
  - Tuple containing the score matrix and the aligned sequences `(aligned_seq1, aligned_seq2)`.

### 2️⃣ Local Alignment (Smith-Waterman)
**File:** `LocalAlignment.py`

This script performs **local sequence alignment** using the Smith-Waterman algorithm, which is useful for finding the most similar subsequences within two sequences.

#### 🔹 Function: `LocalAlignment(seq1, seq2, score)`
- **Parameters:**
  - `seq1` (str): The first sequence.
  - `seq2` (str): The second sequence.
- **Returns:**
  - Tuple containing the aligned subsequences `(aligned_seq1, aligned_seq2)`.

### 3️⃣ BLOSUM Matrix Calculation
**File:** `BLOSUM_calculation.py`

This script calculates the **BLOSUM substitution matrix** based on a given set of sequences, which is used for scoring alignments in sequence comparison.

#### 🔹 Function: `BLOSUM_calculation(sequences)`
- **Parameters:**
  - `sequences` (list of str): A list of sequences to compute the BLOSUM matrix.
- **Returns:**
  - `numpy.array`: The computed BLOSUM matrix.
