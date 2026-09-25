from Bio.Align import substitution_matrices

# global variables
blosum_62 = substitution_matrices.load("BLOSUM62")
gap_penalty = 8
match = 0
go_left = 1
go_up = 2

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
    pointer_matrix = [[-1 for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]

    # Initialising matrices with gap penalties
    for j in range(1, len(seq1) + 1):
        alignment_matrix[0][j] = j * - gap_penalty
        pointer_matrix[0][j] = go_left
    for i in range(1, len(seq2) + 1):
        alignment_matrix[i][0] = i * - gap_penalty
        pointer_matrix[i][0] = go_up

    # Building alignment score and pointer matrices based on Needleman-Wunsch calculations
    for i in range(1, len(seq2) + 1):
        for j in range(1, len(seq1) + 1):
            diagonal = alignment_matrix[i - 1][j - 1] + scoring_function(seq1[j - 1], seq2[i - 1])
            left = alignment_matrix[i][j - 1] - gap_penalty
            up = alignment_matrix[i - 1][j] - gap_penalty
            
            if diagonal >= left and diagonal >= up:
                alignment_matrix[i][j] = diagonal
                pointer_matrix[i][j] = match

            elif left >= up:
                alignment_matrix[i][j] = left
                pointer_matrix[i][j] = go_left

            else:
                alignment_matrix[i][j] = up
                pointer_matrix[i][j] = go_up

    # Traceback and calculating final score
    i = len(seq2)
    j = len(seq1)
    traceback_seq1 = '' 
    traceback_seq2 = ''
    final_score = 0
    
    while (i != 0 or j != 0):
        final_score += alignment_matrix[i][j]

        next = pointer_matrix[i][j]

        if next == match:
            traceback_seq1 += seq1[j - 1]
            traceback_seq2 += seq2[i - 1]
            i = i - 1
            j = j - 1

        elif next == go_left:
            traceback_seq1 += seq1[j - 1]
            traceback_seq2 += '-'
            j = j - 1

        elif next == go_up:
            traceback_seq1 += '-'
            traceback_seq2 += seq2[i - 1]
            i = i - 1
      
            
    aligned_seq1 = traceback_seq1[::-1]
    aligned_seq2 = traceback_seq2[::-1]
         
    return (aligned_seq1, aligned_seq2, final_score)


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
    # 2 matrices to keep track of functions (already initialised)
    alignment_matrix = [[0 for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    pointer_matrix = [[-1 for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]

    max_score = 0
    max_coord = (0, 0)
    # Building alignment score and pointer matrices based on Smith-Waterman calculations
    for i in range(1, len(seq2) + 1):
        for j in range(1, len(seq1) + 1):
            diagonal = alignment_matrix[i - 1][j - 1] + scoring_function(seq1[j - 1], seq2[i - 1])
            left = alignment_matrix[i][j - 1] - gap_penalty
            up = alignment_matrix[i - 1][j] - gap_penalty

            if diagonal <= 0 and left <= 0 and up <= 0:
                continue

            elif diagonal >= left and diagonal >= up:
                alignment_matrix[i][j] = diagonal
                pointer_matrix[i][j] = match

            elif left >= up:
                alignment_matrix[i][j] = left
                pointer_matrix[i][j] = go_left

            else:
                alignment_matrix[i][j] = up
                pointer_matrix[i][j] = go_up

            if alignment_matrix[i][j] >= max_score:
                max_score = alignment_matrix[i][j]
                max_coord = (i, j)
    
    # Traceback and calculating final score
    i, j = max_coord
    traceback_seq1 = '' 
    traceback_seq2 = ''
    final_score = 0
    
    while (alignment_matrix[i][j] != 0):
        final_score += alignment_matrix[i][j]

        next = pointer_matrix[i][j]

        if next == match:
            traceback_seq1 += seq1[j - 1]
            traceback_seq2 += seq2[i - 1]
            i = i - 1
            j = j - 1

        elif next == go_left:
            traceback_seq1 += seq1[j - 1]
            traceback_seq2 += '-'
            j = j - 1

        elif next == go_up:
            traceback_seq1 += '-'
            traceback_seq2 += seq2[i - 1]
            i = i - 1
        
            
    aligned_seq1 = traceback_seq1[::-1]
    aligned_seq2 = traceback_seq2[::-1]
            
    return (aligned_seq1, aligned_seq2, final_score)


## Scoring Function using BLOSUM62
def blosum62_scoring_function(aa_i,aa_j):
    score = blosum_62[aa_i][aa_j]
    return (score)
