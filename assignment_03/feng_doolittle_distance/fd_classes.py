"""
classes used in the code
"""

class pair_alignment:
    """
    class for a pair alignment
    initializes with two seqences and the scoring matrix s
    and gap penalty d
    """

    def __init__(self, seq1, seq2, s, d) -> None:
        self.first = seq1
        self.second = seq2
        self.s = s
        self.d = d


    def compute(self):
        '''
        gets the two sequences seq1 and seq2, 
        the scoring matrix s with match and mismatch scores and
        the gap penalty score d
        computes the dynamic programming matrix with the Needleman-Wunsch algorithm
        returns the dynamic programming matrix
        in this matrix the traceback is included
        '''

        # length of first sequence
        n = len(self.first)

        # length of second sequnce
        m = len(self.second)

        # initialize matrix array
        matrix = []
        for j in range(m+1):
            matrix.append([])
            for i in range(n+1):
                matrix[j].append([0,[0,0]])

        # initialize matrix values
        for i in range(n+1):
                matrix[0][i] = [-i * self.d, [i-1, 0]]
        for j in range(m+1):
            matrix[j][0] = [-j * self.d, [0, j-1]]

        # iterates through all fields of the matrix
        for i in range(1,n+1):
            for j in range(1,m+1):

                # calculate if match or mismatch or gap-gap
                if(self.first[i-1] == "-" and self.second[j-1] == "-"):
                    score = 0
                elif(self.first[i-1] == self.second[j-1]):
                    score = self.s[0]
                else:
                    score = self.s[1]

                # calculate maximum
                maximum = max(matrix[j-1][i-1][0] + score, matrix[j-1][i][0]-self.d, matrix[j][i-1][0]-self.d)

                # sets the max value in the matrix and the traceback
                if(maximum == matrix[j-1][i-1][0] + score):
                    matrix[j][i] = [maximum, [i-1,j-1]]
                elif(maximum == matrix[j-1][i][0]-self.d):
                    matrix[j][i] = [maximum, [i,j-1]]
                else:
                    # maximum == matrix[j][i-1][0]-d could be condition here
                    matrix[j][i] = [maximum, [i-1,j]]

        return matrix

    def getScore(self):
        '''
        returns the score of the global alignment
        '''
        matrix = self.compute()
        return matrix[len(matrix)-1][len(matrix[0])-1][0]

    
    def getAlignment(self):
        '''
        calculates one alignment with the matrix
        '''

        matrix = self.compute()

        # initializing some variables
        start_i = len(matrix)-1
        start_j = len(matrix[0])-1

        seq1_aligned = ""
        seq2_aligned = ""

        old_i = len(matrix)-1
        old_j = len(matrix[0])-1

        count = 0

        # going from (n,m) backwards until we are at (0,0)
        while(start_j != -1 and start_i != -1):
            
            if(old_i == start_i and (count != 0)):
                seq2_aligned += "-"
                seq1_aligned += self.first[old_j-1]
            elif(old_j == start_j and (count != 0)):
                seq1_aligned += "-"
                seq2_aligned += self.second[old_i-1]
            elif(count != 0):
                seq1_aligned += self.first[old_j-1]
                seq2_aligned += self.second[old_i-1]

            count += 1
            old_j = start_j
            old_i = start_i

            start_i = matrix[old_i][old_j][1][1]
            start_j = matrix[old_i][old_j][1][0]

        # reversing the strings
        seq1_aligned = seq1_aligned[::-1]
        seq2_aligned = seq2_aligned[::-1]

        return [seq1_aligned, seq2_aligned]


    def getSobs(self):
        '''
        returns the S_obs Score for the two sequences
        '''
        return self.getScore()

    
    def getSid(self):
        '''
        returns the S_id Score for the two sequences
        '''
        id_seq1 = pair_alignment(self.first, self.first, self.s, self.d)
        id_seq2 = pair_alignment(self.second, self.second, self.s, self.d)

        score_id_seq1 = id_seq1.getScore()
        score_id_seq2 = id_seq2.getScore()

        return (score_id_seq1+score_id_seq2)/2


    def getSrand(self):
        '''
        returns the S_rand Score for the two Sequences
        '''

        # FIXME check if the last part of the formula to calculate the S_rand Score
        # - N(g)*d is inside the sums or is substracted at the end after the sum!!!

        # first sequence
        X = self.first

        # second sequence
        Y = self.second

        # gap penalty
        d = self.d

        # mismatch and match scores
        match = self.s[0]
        mismatch = self.s[1]

        # get the two optimal aligned sequences
        aligned_sequences = self.getAlignment()

        # length of aligned sequences
        L = len(aligned_sequences[0]) * 2

        # calculate number of gaps in optimal alignment
        Ng = 0
        for i in range(len(aligned_sequences[0])):
            if(aligned_sequences[0][i] == "-"):
                Ng += 1
        for i in range(len(aligned_sequences[1])):
            if(aligned_sequences[1][i] == "-"):
                Ng += 1

        # calculate how often every residue appears in the two sequences
        dict_N = {}
        dict_N[X] = {}
        dict_N[Y] = {}

        for i in range(len(X)):
            try:
                dict_N[X][X[i]] += 1
            except:
                dict_N[X][X[i]] = 1
        for i in range(len(Y)):
            try:
                dict_N[Y][Y[i]] += 1
            except:
                dict_N[Y][Y[i]] = 1

        # calculate S_rand
        Srand = 0
        for i in range(len(X)):
            for j in range(len(Y)):
                if(X[i] == Y[j] == "-"):
                    score = 0
                elif(X[i] == Y[j]):
                    score = match
                else:
                    score = mismatch

                Srand += (
                    score *
                    dict_N[X][X[i]] *
                    dict_N[Y][Y[j]]
                )

        Srand = (1/L) * (Srand) - Ng * d

        # FIXME soll Ng*d in jedem summanden subtrahiert werden oder nur einmal zum Schluss?

        return Srand


