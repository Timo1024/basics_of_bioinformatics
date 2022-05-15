
class pair_alignment:
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
        returns array with the new strings (aligned seq1 and seq2)
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

    def getAlignmentPa(self, pa1, pa2):
        '''
        returns array with the new strings (aligned seq1 and seq2)
        and the one new pa which is the combined pa of pa1 and pa2
        with the inserted gaps, which occured when aligning seq1 and seq2
        '''
        matrix = self.compute()

        # save the array of seqs of pa1 and pa2
        seqs1 = pa1.seqs
        seqs2 = pa2.seqs

        # initializing some variables
        start_i = len(matrix)-1
        start_j = len(matrix[0])-1

        seq1_aligned = ""
        seq2_aligned = ""

        seqs1_aligned = []
        seqs2_aligned = []

        for k in seqs1:
            seqs1_aligned.append("")
        for l in seqs2:
            seqs2_aligned.append("")

        old_i = len(matrix)-1
        old_j = len(matrix[0])-1

        count = 0

        # going from (n,m) backwards until we are at (0,0)
        while(start_j != -1 and start_i != -1):
            
            if(old_i == start_i and (count != 0)):
                seq2_aligned += "-"
                for i in range(len(seqs2)):
                    seqs2_aligned[i] += "-"

                seq1_aligned += self.first[old_j-1]
                for i in  range(len(seqs1)):
                    seqs1_aligned[i] += seqs1[i][old_j-1]

            elif(old_j == start_j and (count != 0)):
                seq1_aligned += "-"
                for i in range(len(seqs1)):
                    seqs1_aligned[i] += "-"

                seq2_aligned += self.second[old_i-1]
                for i in  range(len(seqs2)):
                    seqs2_aligned[i] += seqs2[i][old_i-1]

            elif(count != 0):
                seq1_aligned += self.first[old_j-1]
                for i in  range(len(seqs1)):
                    seqs1_aligned[i] += seqs1[i][old_j-1]

                seq2_aligned += self.second[old_i-1]
                for i in  range(len(seqs2)):
                    seqs2_aligned[i] += seqs2[i][old_i-1]

            count += 1
            old_j = start_j
            old_i = start_i

            start_i = matrix[old_i][old_j][1][1]
            start_j = matrix[old_i][old_j][1][0]

        # reversing the strings
        seq1_aligned = seq1_aligned[::-1]
        seq2_aligned = seq2_aligned[::-1]

        # reversing all strings in the two seqs arrays
        for i in range(len(seqs1_aligned)):
            seqs1_aligned[i] = seqs1_aligned[i][::-1]
        for i in range(len(seqs2_aligned)):
            seqs2_aligned[i] = seqs2_aligned[i][::-1]

        # combining the two arrays in one
        seqs_combined = []
        for i in seqs1_aligned:
            seqs_combined.append(i)
        for i in seqs2_aligned:
            seqs_combined.append(i)

        return [seq1_aligned, seq2_aligned, seqs_combined]


class profile_alignment:
    def __init__(self, aligned_seqs) -> None:
        self.seqs = aligned_seqs

    def getID(self):
        # combines all sequences so that in calc_scross 
        # the used pas can be found and be removed
        return "".join(self.seqs)

    def getGapped(self, seq):
        '''
        returns the aligned_seqs with the inserted gap, where seq has the gaps
        '''
        # TODO use old seq to see where the new gaps are b/c there could be gaps before the alignment too
        # init new array
        seqs_gapped = []
        for x in self.seqs:
            seqs_gapped.append("")
        for i in range(len(seq)):
            if(seq[i] != "-"):
                for k in range(len(seqs_gapped)):
                    seqs_gapped[k] += seq[i]
            else:
                for k in range(len(seqs_gapped)):
                    seqs_gapped[k] += "-"
        return seqs_gapped