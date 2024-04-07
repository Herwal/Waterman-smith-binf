# BINF100 2020
# code written by Susanna Roeblitz
# This is an exmample code that realizes the Smith-Waterman algorithm for local alignment.
# Note that the alignment found is very sensitive to the chosen scoring scheme.
# For example, changing the gap penalty from 2 to 1 changes drastically the
# number of optimal local alignments in the virus example.

import numpy as np

# Step 1: read in sequences
seq1 = ""
with open("caliSeq2009.fasta", "r") as infile1:
    # ignore FASTA comment
    infile1.readline()
    lines = infile1.readlines(100000)  # read 100,000 bytes
    while lines:
        seq1 = seq1 + "".join(lines)
        lines = infile1.readlines(100000)
    # replace newline character in joined string
    seq1 = seq1.replace("\n", "")


seq2 = ""
with open("brisbaneSeq2007.fasta", "r") as infile2:
    # ignore FASTA comment
    infile2.readline()
    lines = infile2.readlines(100000)  # read 100,000 bytes
    while lines:
        seq2 = seq2 + "".join(lines)
        lines = infile2.readlines(100000)
    # replace newline character in joined string
    seq2 = seq2.replace("\n", "")

# test sequences:
#
# seq1='AAAGCTCCGATCTCG'
# seq2='TAAAGCAATTTTGGTTTTTTTCCGA'

# seq1='AGGGCT'
# seq2='AGGCA'

# seq1='AAACTGTCAAA'
# seq2='AAACTTCAAA'

# seq1='AAACTGTCAAA'
# seq2='AAACTCAAA'

# seq1='AAACTTCAAA'
# seq2='AAACTCAAA'

# Initialization
N = len(seq1)
M = len(seq2)
# initialize (M+1) by (N+1) matrix
H = np.zeros((M + 1, N + 1))

# scoring scheme
gap = -2
mismatch = -1
match = 1

# Step 1: Build Alignment matrix
# H[0][0]=0
for j in range(1, N + 1, 1):
    # initialize 1st row with zeros (free initial gaps in string 2)
    H[0][j] = 0
for i in range(1, M + 1, 1):
    # initialize 1st column with zeros (free initial gaps in string 1)
    H[i][0] = 0
for i in range(1, M + 1, 1):
    for j in range(1, N + 1, 1):
        if seq1[j - 1] == seq2[i - 1]:
            score1 = H[i - 1][j - 1] + match
        else:
            score1 = H[i - 1][j - 1] + mismatch
        score2 = H[i][j - 1] + gap
        score3 = H[i - 1][j] + gap
        # take max of candidate scores and zero:
        H[i][j] = max(0, score1, score2, score3)

# print(H)


def printAlignment(i, j, A1, A2):
    matches = 0
    matchstr = ""
    if len(A1) != len(A2):
        print("Something is wrong with the alignment length!")
    lenAlignment = len(A1)
    # count number of matches:
    for s in range(0, lenAlignment, 1):
        if A1[s] == A2[s]:
            matches = matches + 1
            matchstr = matchstr + "|"
        else:
            matchstr = matchstr + " "
    # print('sequence identity:',round(matches/min(len(seq1),len(seq2)),2))
    print("sequence identity:", 100 * round(matches / lenAlignment, 3), "%")
    a1 = A1[::-1]
    a2 = A2[::-1]
    matchstr = matchstr[::-1]
    # split lines into blocks of length blocklen
    blocklen = 50
    blocks = int(lenAlignment / blocklen) + 1
    posA1 = j + 1
    posA2 = i + 1
    for bl in range(1, blocks, 1):
        # extract substrings:
        subA1 = a1[((bl - 1) * blocklen) : (bl * blocklen)]
        subA2 = a2[((bl - 1) * blocklen) : (bl * blocklen)]
        submatch = matchstr[((bl - 1) * blocklen) : (bl * blocklen)]
        # count number of gaps:
        ngA1 = len(subA1.split("-")) - 1
        ngA2 = len(subA2.split("-")) - 1
        # print lines:
        print("seq1", posA1, subA1, posA1 + blocklen - 1 - ngA1, sep="\t")
        print("\t", " " * 6, submatch)
        print("seq2", posA2, subA2, posA2 + blocklen - 1 - ngA2, sep="\t")
        print("\n")
        posA1 = posA1 + blocklen - ngA1
        posA2 = posA2 + blocklen - ngA2
    # treat last block special since it is shorter
    subA1 = a1[((blocks - 1) * blocklen) : lenAlignment]
    subA2 = a2[((blocks - 1) * blocklen) : lenAlignment]
    submatch = matchstr[((blocks - 1) * blocklen) : lenAlignment]
    L = len(subA1)
    # count number of gaps:
    ngA1 = len(subA1.split("-")) - 1
    ngA2 = len(subA2.split("-")) - 1
    # print lines:
    print("seq1", posA1, subA1, posA1 + L - ngA1 - 1, sep="\t")
    print("\t", " " * 6, submatch)
    print("seq2", posA2, subA2, posA2 + L - ngA2 - 1, sep="\t")


def backTrack(i, j, A1, A2, k):
    # i: row index
    # j: column index
    if H[i][j] == 0:
        print("found optimal alignment:")
        printAlignment(i, j, A1, A2)
    else:
        # vertical move
        if H[i][j] == (H[i - 1][j] + gap):
            A1mod = A1 + "-"
            A2mod = A2 + seq2[i - 1]
            backTrack(i - 1, j, A1mod, A2mod, k + 1)
        # horizontal move
        if H[i][j] == (H[i][j - 1] + gap):
            A1mod = A1 + seq1[j - 1]
            A2mod = A2 + "-"
            backTrack(i, j - 1, A1mod, A2mod, k + 1)
        # diagonal move
        if (
            (H[i][j] == (H[i - 1][j - 1] + match)) and (seq2[i - 1] == seq1[j - 1])
        ) or (
            (H[i][j] == (H[i - 1][j - 1] + mismatch)) and (seq2[i - 1] != seq1[j - 1])
        ):
            A1mod = A1 + seq1[j - 1]
            A2mod = A2 + seq2[i - 1]
            backTrack(i - 1, j - 1, A1mod, A2mod, k + 1)


# Find the largest value(s) and its(their) index in the matrix
result = np.where(H == np.amax(H))

# print('Tuple of arrays returned : ', result)

# zip the 2 arrays to get the exact coordinates
listOfCordinates = list(zip(result[0], result[1]))

print("optimal local alignment score:", np.amax(H))

# backtrack from all starting positions
k = 0
# traverse over the list of cordinates
for cord in listOfCordinates:
    print("start backtracking at position:")
    print(cord)

    A1 = ""
    A2 = ""

    rowIdx = int(cord[0])
    colIdx = int(cord[1])

    backTrack(rowIdx, colIdx, A1, A2, len(A1))
    k = k + 1
