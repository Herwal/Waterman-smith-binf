from pathlib import Path

# THIS ONE WORKS !!!

# constants
GAP_PENALTY = -1
MATCH = 1
MISMATCH = 0

# local alignment
INITIAL_GAP_PENALTY = 0

# num of align with highest score
ALIGNMENT_COUNTER = 0


def create_table(SEQ1: str, SEQ2: str) -> list[list[int]]:
    """
    Creates a dynamic programming table for sequence alignment.

    Args:
        SEQ1 (str): The first sequence.
        SEQ2 (str): The second sequence.

    Returns:
        list: A 2D list representing the dynamic programming table.
    """

    # empty table
    table = [[None for _ in range(len(SEQ1) + 1)] for _ in range(len(SEQ2) + 1)]

    # always 0
    table[0][0] = 0

    # setup minus for first row
    for i in range(len(table[0])):
        table[0][i] = INITIAL_GAP_PENALTY * i

    # setup minus for all rows
    for i in range(len(table)):
        table[i][0] = GAP_PENALTY * i

    # loop through table
    for row_index in range(1, len(table)):
        for col_index in range(1, len(table[0])):

            seq1_base = SEQ1[col_index - 1]
            seq2_base = SEQ2[row_index - 1]

            # check above
            ABOVE_SCORE = table[row_index - 1][col_index] + GAP_PENALTY

            # check left
            LEFT_SCORE = table[row_index][col_index - 1] + GAP_PENALTY

            # check diag
            DIAG_SCORE = table[row_index - 1][col_index - 1]
            # if dna base match
            if seq1_base == seq2_base:
                DIAG_SCORE += MATCH
            else:
                DIAG_SCORE += MISMATCH

            table[row_index][col_index] = max(ABOVE_SCORE, LEFT_SCORE, DIAG_SCORE)

    # print(table)
    # for i in table:
    #     print(i)

    return table


def alignment_score(align1: str, align2: str) -> int:
    """
    Calculates the score of an alignment.

    Args:
        align1 (str): The first aligned sequence.
        align2 (str): The second aligned sequence.

    Returns:
        int: The score of the alignment.
    """
    score = 0
    for a, b in zip(align1, align2):
        if a == "-" or b == "-":
            score += GAP_PENALTY
        elif a == b:
            score += MATCH
        else:
            score += MISMATCH
    return score


def find_all_possible_alignments(
    i: int, j: int, A1: str, A2: str, table: list[list], max_score
):
    """
    Recursively finds all possible alignments between two sequences based on a scoring table.

    Args:
        i (int): The current row index in the scoring table.
        j (int): The current column index in the scoring table.
        A1 (str): The partial alignment of sequence 1.
        A2 (str): The partial alignment of sequence 2.
        table (list[list]): The scoring table.
        max_score: The maximum alignment score.

    Returns:
        tuple: A tuple containing the reversed alignments of sequence 1 and sequence 2.

    """
    if (i > 0) and (j > 0):
        match_score = MATCH if SEQ1[j - 1] == SEQ2[i - 1] else MISMATCH

        # check above
        if i > 0 and table[i][j] == table[i - 1][j] + GAP_PENALTY:
            find_all_possible_alignments(
                i - 1, j, A1 + "-", A2 + SEQ2[i - 1], table, max_score
            )

        # check left
        if j > 0 and table[i][j] == table[i][j - 1] + GAP_PENALTY:
            find_all_possible_alignments(
                i, j - 1, A1 + SEQ1[j - 1], A2 + "-", table, max_score
            )

        # check diag
        if i > 0 and j > 0 and table[i][j] == table[i - 1][j - 1] + match_score:
            find_all_possible_alignments(
                i - 1, j - 1, A1 + SEQ1[j - 1], A2 + SEQ2[i - 1], table, max_score
            )

    else:
        alignment_score_current = alignment_score(A1[::-1], A2[::-1])
        if alignment_score_current >= max_score:
            print(
                f"Alignment: Score: {alignment_score_current} {A1[::-1]} | {A2[::-1]}"
            )
            global ALIGNMENT_COUNTER
            ALIGNMENT_COUNTER += 1

        return A1[::-1], A2[::-1]


def get_max_score_indices(table: list[list[int]]) -> tuple[list[tuple[int, int]], int]:
    """
    Returns the indices of the maximum score(s) in the given table.

    Args:
        table (list): A 2D list representing the table of scores.

    Returns:
        tuple: A tuple containing a list of indices and the maximum score.
    """
    max_score = float("-inf")
    max_indices = []
    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] == max_score:
                max_indices.append((i, j))
            elif table[i][j] > max_score:
                max_score = table[i][j]
                max_indices = [(i, j)]
    return max_indices, max_score


def find_best_alignment(i: int, j: int, A1: str, A2: str, table: list[list]) -> int:
    """
    Recursively finds the best alignment score between two sequences.

    Args:
        i (int): The current row index in the table.
        j (int): The current column index in the table.
        A1 (str): The partial alignment of sequence 1.
        A2 (str): The partial alignment of sequence 2.
        table (list[list]): The scoring table.

    Returns:
        int: The best alignment score.

    """
    if (i > 0) and (j > 0):
        print("finding alingment...")
        match_score = MATCH if SEQ1[j - 1] == SEQ2[i - 1] else MISMATCH

        # check above
        if i > 0 and table[i][j] == table[i - 1][j] + GAP_PENALTY:
            return find_best_alignment(i - 1, j, A1 + "-", A2 + SEQ2[i - 1], table)

        # check left
        if j > 0 and table[i][j] == table[i][j - 1] + GAP_PENALTY:
            return find_best_alignment(i, j - 1, A1 + SEQ1[j - 1], A2 + "-", table)

        # check diag
        if i > 0 and j > 0 and table[i][j] == table[i - 1][j - 1] + match_score:
            return find_best_alignment(
                i - 1, j - 1, A1 + SEQ1[j - 1], A2 + SEQ2[i - 1], table
            )
    else:
        score = alignment_score(A1[::-1], A2[::-1])
        print(f"BEST Alignment: Score: {score} {A1[::-1]} | {A2[::-1]}")
        return score


if __name__ == "__main__":

    SEQ1 = "AAAGCTCCGATCTCG"
    SEQ2 = "TAAAGCAATTTTGGTTTTTTTCCGA"

    # SEQ1 = Path("resources/DNA_A_California_2009_pandemicH1N1_segment7.txt").read_text().replace('\n', '')
    # SEQ2 = Path("resources/DNA_A_Brisbane_2007_H1N1_M2_CDS.txt").read_text().replace('\n', '')

    table = create_table(SEQ1, SEQ2)

    A1 = ""
    A2 = ""
    max_score_indices, max_score = get_max_score_indices(table)
    best_alignment_score = find_best_alignment(
        max_score_indices[0][0], max_score_indices[0][1], A1, A2, table
    )
    for max_i, max_j in max_score_indices:
        find_all_possible_alignments(max_i, max_j, A1, A2, table, best_alignment_score)

    # look through entire table
    # for i in range(len(table)):
    #     for j in range(len(table[0])):
    #         find_all_possible_alignments(i, j, A1, A2, table, best_alignment_score)

    print(f"Amount of alignments found with score {max_score} = {ALIGNMENT_COUNTER}")
