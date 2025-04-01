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

The **Needleman-Wunsch algorithm** is a dynamic programming method used for **global sequence alignment**. It ensures that two sequences are aligned from start to finish by introducing gaps where necessary. This algorithm is widely used in genomics and proteomics for comparing DNA, RNA, or protein sequences.

🔹 **How It Works:**
- A scoring matrix is initialized based on match, mismatch, and gap penalties.
- The matrix is filled using dynamic programming to find the optimal alignment.
- A traceback step reconstructs the best alignment path.

🔹 **Applications:**
- Comparing full-length gene sequences.
- Identifying evolutionary relationships between species.
- Aligning protein sequences for structural analysis.

---

### 2️⃣ Local Alignment (Smith-Waterman)
**File:** `LocalAlignment.py`

The **Smith-Waterman algorithm** is a dynamic programming approach used for **local sequence alignment**. Unlike global alignment, this method finds **only the most similar subsections** of two sequences, making it ideal for identifying conserved regions within genes or proteins.

🔹 **How It Works:**
- A scoring matrix is built, but scores below zero are set to zero (allowing local matches to stand out).
- The highest-scoring region is identified for optimal local alignment.
- A traceback step extracts the aligned subsequences.

🔹 **Applications:**
- Detecting functional domains in proteins.
- Finding gene sequences within larger genomic data.
- Comparing short DNA or protein segments to reference databases.

---

### 3️⃣ BLOSUM Matrix Calculation
**File:** `BLOSUM_calculation.py`

The **BLOSUM (BLOcks Substitution Matrix)** is a scoring matrix used in sequence alignment for comparing protein sequences. It is derived from aligned blocks of protein sequences that share evolutionary relationships.

🔹 **How It Works:**
- A set of related protein sequences is grouped into clusters based on sequence identity.
- The frequency of amino acid substitutions is calculated to generate substitution scores.
- The resulting matrix helps determine the likelihood of one amino acid mutating into another.

🔹 **Applications:**
- Used in BLAST (Basic Local Alignment Search Tool) for protein sequence comparison.
- Helps in predicting protein function and evolutionary distance.
- Aids in identifying conserved residues important for protein stability.
