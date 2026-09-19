from Bio import substitution_matrices
def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    # 2 matrices to keep track of functions
    alignment_matrix = [[0 for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    pointer_matrix = [[None for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]

    # Initialising matrices with gap penalties (assuming 1)
    for j in range(1, len(seq1) + 1):
        alignment_matrix[0][j] = j * - 1
        pointer_matrix[0][j] = (0, j - 1)
    for i in range(1, len(seq2) + 1):
        alignment_matrix[i][0] = i * -1
        pointer_matrix[i][0] = (i - 1, 0)

    # Building alignment score and pointer matrices based on Needleman-Wunsch calculations
    for i in range(2, len(seq2) + 1):
        for j in range(2, len(seq1) + 1):
            alignment_matrix[i][j] = max(
                (alignment_matrix[i - 1][j - 1] + scoring_function(seq1[j - 1], seq2[i - 1])), 
                (alignment_matrix[i - 1][j] - 1),
                (alignment_matrix[i][j - 1] - 1))
            
            if alignment_matrix[i][j] == (alignment_matrix[i - 1][j - 1] + scoring_function(seq1[j - 1], seq2[i - 1])):
                pointer_matrix[i][j] = (i - 1, j - 1)
            elif alignment_matrix[i][j] == (alignment_matrix[i - 1][j] - 1):
                pointer_matrix[i][j] = (i - 1, j)
            elif alignment_matrix[i][j] == (alignment_matrix[i][j - 1] - 1):
                pointer_matrix[i][j] = (i, j - 1)

    # NEED TO IMPLEMENT: 
    # - Traceback + building final aligned sequences
    # - Final alignment score

         
    raise NotImplementedError()


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    raise NotImplementedError()


## Scoring Function using BLOSUM62
def scoring_function(aa_i,aa_j):
    blosum_62 = substitution_matrices.load("BLOSUM62")
    score = blosum_62[aa_i][aa_j]
    return (score)
