"""
Implementasjon av Smith-Waterman algoritmen for lokal sammenligning av to sekvenser. WIP.
"""

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
    """
    Function to read a fasta file and return the sequence as a string.
    Skips the first line of the file.
    @param file_name: The name of the file to read.
    """
    with open(file_name, "r", encoding="utf-8") as f:
        f.readline()
        return f.read().replace("\n", "")


def print_alignment(seq1: str, seq2: str, start_positions: list) -> None:
    """
    Function to print the alignment of two sequences in a pretty format.
    @constant MAX_WIDTH: The maximum width of the alignment, length of each line.
    @param seq1: The first sequence.
    @param seq2: The second sequence.
    """
    MAX_WIDTH = 50
    start_row = start_positions[0] - len(seq1)
    start_col = start_positions[1] - len(seq2)

    for i in range(0, len(seq1), MAX_WIDTH):
        # Calculate the number of spaces to print before the sequence
        num_gap_upper = " " * (len(str(start_col)) - 3)
        num_gap_lower = " " * (len(str(start_row)) - 3)
        # Whitespace gap is just the length of the word "Sequence1: " + the length of the number of the start row
        white_space_gap = "            " + " " * (len(str(start_row)))

        # Sequence 1
        print(
            f"Sequence1: {start_row} {num_gap_upper}{seq1[i:i+MAX_WIDTH]} {start_row + len(seq1[i:i+MAX_WIDTH])}"
        )
        # '|' if match, ' ' if mismatch
        print(
            white_space_gap
            + "".join(
                "|" if seq1[i + k] == seq2[i + k] else " "
                for k in range(min(MAX_WIDTH, len(seq1) - i))
            )
        )
        # Sequence 2
        print(
            f"Sequence2: {start_col} {num_gap_lower}{seq2[i:i+MAX_WIDTH]} {start_col + len(seq2[i:i+MAX_WIDTH])}"
        )
        print()
        start_row += MAX_WIDTH
        start_col += MAX_WIDTH


def create_alignment_matrix(seq1: str, seq2: str) -> np.ndarray:
    """
    Function to create the alignment matrix for two given sequences.
    Uses numpy.zeros to initialize the matrix with m + 1 and n + 1 dimensions.
    @param seq1: The first sequence.
    @param seq2: The second sequence.
    @return: The alignment matrix.
    """
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
    """
    Function to perform local alignment of two sequences.
    Uses create_alignment_matrix and backtrack functions to find the optimal alignment(s).
    Uses print_statistics to print the statistics for each optimal alignment.
    @param seq1: The first sequence.
    @param seq2: The second sequence.
    """

    # Create the alignment matrix and get the optimal alignment score
    aligned_sequences = []
    start_positions = []
    matrix = create_alignment_matrix(seq1, seq2)
    alignment_score = np.max(matrix)
    max_pos = np.where(matrix == alignment_score)

    # Backtrack to find the optimal alignment(s) for each alignment tuple
    for start_row, start_col in zip(*max_pos):
        # Calles the params class to hold the parameters for the backtrack function
        back_track_params = BacktrackParams(
            start_row, start_col, "", "", aligned_sequences, matrix, seq1, seq2
        )
        backtrack(back_track_params)
        start_positions.append((start_row, start_col))

    print("Matrix: \n", matrix, "\n")
    print("Optimal alignment score:", alignment_score)
    print(f"There are {len(aligned_sequences)} optimal alignment(s).")

    # Get and print the statistics for each optimal alignment
    for idx, (a1, a2) in enumerate(aligned_sequences):
        length = len(a1)
        # Get the number of matches, mismatches and gaps in the alignment
        matches, mismatches, gaps = get_alignment_statistics(aligned_sequences[idx])

        # Calls the params class to hold the parameters for the print statistics function
        stats_params = StatisticsParams(
            matches, mismatches, gaps, length, aligned_sequences[idx], idx
        )
        # Calls print statistics function
        print_statistics(stats_params, tuple(start_positions[idx]))


def backtrack(params: BacktrackParams) -> None:
    """
    Function to backtrack and find the optimal alignment(s) using recursion.
    Uses the defined scoring scheme to determine the best move.
    @param params: The parameters for the backtrack function.
    @return: The optimal alignment(s).
    """
    row = params.row
    col = params.col
    a1 = params.a1
    a2 = params.a2
    alignments = params.alignments
    matrix = params.matrix
    seq1 = params.seq1
    seq2 = params.seq2

    # Base case
    if matrix[row][col] == 0:
        alignments.append((a1[::-1], a2[::-1]))
        return alignments

    # Check if it's a horizontal move
    if col > 0 and matrix[row][col] == (matrix[row][col - 1] + GAP):
        new_params = BacktrackParams(
            row,
            col - 1,
            a1 + seq1[col - 1],
            a2 + "-",
            alignments,
            matrix,
            seq1,
            seq2,
        )
        return backtrack(new_params)

    # Check if it's a vertical move
    if row > 0 and matrix[row][col] == (matrix[row - 1][col] + GAP):
        new_params = BacktrackParams(
            row - 1,
            col,
            a1 + "-",
            a2 + seq2[row - 1],
            alignments,
            matrix,
            seq1,
            seq2,
        )
        return backtrack(new_params)

    # Check if it's a diagonal move (match/mismatch)
    if (
        row > 0
        and col > 0
        and matrix[row][col]
        == matrix[row - 1][col - 1]
        + (MATCH if seq1[col - 1] == seq2[row - 1] else MISMATCH)
    ):
        new_params = BacktrackParams(
            row - 1,
            col - 1,
            a1 + seq1[col - 1],
            a2 + seq2[row - 1],
            alignments,
            matrix,
            seq1,
            seq2,
        )
        return backtrack(new_params)


def get_alignment_statistics(aligned_sequences: list) -> tuple:
    """
    Function to calculate the statistics for an alignment.
    @param aligned_sequences: a list of the aligned sequences, contains only 1 alignment for each call.
    @return: The number of matches, mismatches and gaps in the alignment.
    """
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


def print_statistics(params: StatisticsParams, start_positions) -> None:
    """
    Function to print the statistics for an alignment.
    @param params: The parameters for the print statistics function.
    @param start_positions: The start positions of the alignment.
    """
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

    print_alignment(a1, a2, start_positions)


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
