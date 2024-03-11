"""Implementasjon av Smith-Waterman algoritmen for lokal sammenligning av to sekvenser."""

import numpy as np

# Define scoring scheme
MATCH = 1
MISMATCH = -1
GAP = -2
MINIMUM_SCORE = 0

brisbane_seq_2007 = "brisbaneSeq2007.fasta"
cali_seq_2009 = "caliSeq2009.fasta"


class BacktrackParams:
    """Class to hold the parameters for the backtrack function."""

    def __init__(self, row, col, a1, a2, alignments, matrix, seq1, seq2):
        self.row = row
        self.col = col
        self.a1 = a1
        self.a2 = a2
        self.alignments = alignments
        self.matrix = matrix
        self.seq1 = seq1
        self.seq2 = seq2


class StatisticsParams:
    """Class to hold the parameters for the print statistics function."""

    def __init__(self, matches, mismatches, gaps, length, aligned_sequences, idx):
        self.matches = matches
        self.mismatches = mismatches
        self.gaps = gaps
        self.length = length
        self.aligned_sequences = aligned_sequences
        self.idx = idx


def read_fasta(file_name: str) -> str:
    """Function to read a fasta file and return the sequence as a string."""
    with open(file_name, "r", encoding="utf-8") as f:
        f.readline()
        return f.read().replace("\n", "")


def print_alignment(a1: str, a2: str) -> None:
    """ "Function to print the alignment of two sequences."""
    max_width = 20
    for i in range(0, len(a1), max_width):
        end_pos = i + max_width if len(a1) > i + max_width else len(a1)
        print(f"Sequence1: {i + 1: >4} {a1[i:i+max_width]} {end_pos}")
        print(
            "                "
            + "".join(
                "|" if a1[i + k] == a2[i + k] else " "
                for k in range(min(max_width, len(a1) - i))
            )
        )
        print(f"Sequence2: {i + 1: >4} {a2[i:i+max_width]} {end_pos}")
        print()


def create_alignment_matrix(seq1: str, seq2: str) -> np.ndarray:
    """Function to create the alignment matrix for two given sequences."""
    # Calculate the necessary dimensions and initialize the matrix
    N = len(seq1)
    M = len(seq2)
    matrix = np.zeros((M + 1, N + 1))

    # Main loop for the table
    for i in range(1, M + 1):  # rows
        for j in range(1, N + 1):  # cols
            if seq1[j - 1] == seq2[i - 1]:
                diagonal = MATCH + matrix[i - 1][j - 1]

            else:
                diagonal = MISMATCH + matrix[i - 1][j - 1]

            vertical = matrix[i - 1][j] + GAP  # moves one down from previous row
            horisontal = (
                matrix[i][j - 1] + GAP
            )  # moves one to the right from previous col
            matrix[i][j] = max(
                MINIMUM_SCORE, diagonal, vertical, horisontal
            )  # best score for that cell

    return matrix


def local_alignment(seq1: str, seq2: str) -> None:
    """Function to perform local alignment of two sequences."""

    # Create the alignment matrix and get the optimal alignment score
    matrix = create_alignment_matrix(seq1, seq2)
    alignment_score = np.max(matrix)
    max_pos = np.where(matrix == alignment_score)
    start_row, start_col = max_pos[0][0], max_pos[1][0]

    # Start backtracking
    aligned_sequences = []
    back_track_params = BacktrackParams(
        start_row, start_col, "", "", aligned_sequences, matrix, seq1, seq2
    )
    backtrack(back_track_params)

    print("Matrix: \n", matrix, "\n")
    print("Optimal alignment score:", alignment_score)
    print(f"There are {len(aligned_sequences)} optimal alignment(s).")

    # Get and print the statistics for each optimal alignment
    for idx, (a1, a2) in enumerate(aligned_sequences):
        length = len(a1)
        matches, mismatches, gaps = get_alignment_statistics(aligned_sequences[idx])

        stats_params = StatisticsParams(
            matches, mismatches, gaps, length, aligned_sequences[idx], idx
        )
        print_statistics(stats_params)


def backtrack(params: BacktrackParams) -> None:
    """Function to backtrack and find the optimal alignment(s)."""

    # Base case
    if params.matrix[params.row][params.col] == 0:
        params.alignments.append((params.a1[::-1], params.a2[::-1]))

    # Check if it's a horizontal move
    if params.col > 0 and params.matrix[params.row][params.col] == (
        params.matrix[params.row][params.col - 1] + GAP
    ):
        new_params = BacktrackParams(
            params.row,
            params.col - 1,
            params.a1 + params.seq1[params.col - 1],
            params.a2 + "-",
            params.alignments,
            params.matrix,
            params.seq1,
            params.seq2,
        )
        backtrack(new_params)

    # Check if it's a vertical move
    if params.row > 0 and params.matrix[params.row][params.col] == (
        params.matrix[params.row - 1][params.col] + GAP
    ):
        new_params = BacktrackParams(
            params.row - 1,
            params.col,
            params.a1 + "-",
            params.a2 + params.seq2[params.row - 1],
            params.alignments,
            params.matrix,
            params.seq1,
            params.seq2,
        )
        backtrack(new_params)

    # Check if it's a diagonal move (match/mismatch)
    if (
        params.row > 0
        and params.col > 0
        and params.matrix[params.row][params.col]
        == params.matrix[params.row - 1][params.col - 1]
        + (
            MATCH
            if params.seq1[params.col - 1] == params.seq2[params.row - 1]
            else MISMATCH
        )
    ):
        new_params = BacktrackParams(
            params.row - 1,
            params.col - 1,
            params.a1 + params.seq1[params.col - 1],
            params.a2 + params.seq2[params.row - 1],
            params.alignments,
            params.matrix,
            params.seq1,
            params.seq2,
        )
        backtrack(new_params)


def get_alignment_statistics(aligned_sequences: list) -> tuple:
    matches, mismatches, gaps = 0, 0, 0
    a1, a2 = aligned_sequences

    for i in range(len(a1)):
        if a1[i] != "-" and a2[i] != "-":
            if a1[i] == a2[i]:
                matches += 1
            else:
                mismatches += 1
        else:
            gaps += 1

    return matches, mismatches, gaps


def print_statistics(params: StatisticsParams) -> None:
    matches, mismatches, gaps, length, aligned_sequences, idx = (
        params.matches,
        params.mismatches,
        params.gaps,
        params.length,
        params.aligned_sequences,
        params.idx,
    )

    sequence_identity = matches / length

    # Print the optimal alignment
    print(f"\nOptimal Alignment {idx + 1}:")
    print(
        f"Sequence identity: {matches}/{length} ({sequence_identity:.2%}) Mismatches: {mismatches}/{length} ({(mismatches / length):.2%}) Gaps: {gaps}/{length} ({(gaps / length):.2%})\n"
    )
    a1, a2 = aligned_sequences
    print_alignment(a1, a2)


# a)
# test sequences
seq1 = "AAAGCTCCGATCTCG"
seq2 = "TAAAGCAATTTTGGTTTTTTTCCGA"
print("a)\n")
local_alignment(seq1, seq2)


# b)
# File sequences
seq3 = read_fasta(brisbane_seq_2007)
seq4 = read_fasta(cali_seq_2009)
print("\nb)\n")
local_alignment(seq3, seq4)
