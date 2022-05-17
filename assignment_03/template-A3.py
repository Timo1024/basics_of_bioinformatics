#!/usr/bin/env python
"""
Taken from: 

https://gist.github.com/slowkow/06c6dba9180d013dfd82bec217d22eb5

The Needleman-Wunsch Algorithm
==============================
This is a dynamic programming algorithm for finding the optimal alignment of
two strings.
Example
-------
    >>> x = "GATTACA"
    >>> y = "GCATGCU"
    >>> print(nw(x, y))
    G-ATTACA
    GCA-TGCU
LICENSE
This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.
In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
For more information, please refer to <http://unlicense.org/>
"""

from asyncio.windows_events import NULL
import numpy as np
import argparse


def create_parser():
    '''
    Creates a argument parser to facilitate file reading.

    '''
    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-f1', '--file-one',
                   help="File for the the programming task")

    return(p.parse_args())


class alignment:
    def __init__(self, first_seq, second_seq, s, d) -> None:
        self.first = first_seq
        self.second = second_seq
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

                # calculate if match or mismatch
                if(self.first[i-1] == self.second[j-1]):
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
        matrix = self.compute()
        return matrix[len(matrix)-1][len(matrix[0])-1][0]

    def getAlignment(self):
        pass


def extract_headings_sequences(file):
    '''
    saves all headings and sequences in two arrays
    '''

    heading_array  = []
    sequence_array = []
    with open(file) as file_content:
        for line in file_content:
            normal_line = line
            line = line[:line.rfind("\n")]
            if(not normal_line.endswith("\n")):
                line = normal_line
            if(line.startswith(">")):
                heading_array.append(line)
            else:
                if(len(sequence_array) != len(heading_array)):
                    sequence_array.append([])
                sequence_array[len(heading_array)-1].append(line)

    # concats all line which represent one sequence
    sequence_array_concat = []
    for sequence in sequence_array:
        seq = ""
        for line in sequence:
            seq += line
        sequence_array_concat.append(seq)

    return [heading_array, sequence_array_concat]


def find_max_and_rest(all_alignments):
    '''
    finds the alignment with the maximal score from the given alignments in all_alignments
    and gives the alignment of the two other sequences back
    '''

    # get all scores in one array
    scores = []
    for align in all_alignments:
        scores.append(align.getScore())

    # get alignment with max score (A*max)
    max_score_index = np.argmax(scores)
    max_alignment = all_alignments[max_score_index]

    # get alignment of the other two sequences (A*rest)
    a_1 = set([max_alignment.first, max_alignment.second])
    rest_alignment = NULL
    for align in all_alignments:
        a_2 = set([align.first, align.second])
        if not a_1.intersection(a_2):
            rest_alignment = align

    return [max_alignment, rest_alignment]


def find_cross(max_rest):
    '''
    find the A*cross alignment
    '''

    s = max_rest[0][0].s
    d = max_rest[0][0].d

    # make all possible alignments between
    new_alignments = []
    
    pass


def read_arguments():
    '''
    reads the match/mismatch and gap penalty from the command line
    '''

    input_1 = False
    input_2 = False
    input_3 = False

    match = 0
    mismatch = 0
    gap = 0

    # gets input for matching score
    while(not input_1):
        match = input("Type the score for a matching pair:")
        try:
            match = int(match)
            input_1 = True
        except:
            print("Wrong input")

    # gets input for mismatching score
    while(not input_2):
        mismatch = input("Type the score for a mismatch:")
        try:
            mismatch = int(mismatch)
            input_2 = True
        except:
            print("Wrong input")

    # gets input for gap penalty
    while(not input_3):
        gap = input("Type the score for gap penalty:")
        try:
            gap = int(gap)
            input_3 = True
        except:
            print("Wrong input")

    return [match, mismatch, gap]


def main():
    """
    The main function should contain all functions that solved the stated tasks.
    """

    # T2

    # extracting the sequences of the file into one array
    file1 = args.file_one
    sequences = extract_headings_sequences(file1)[1]

    # reads parameters for the alignment
    parameters = read_arguments()

    # scoring matrix [match score, mismatch score]
    s = [parameters[0], parameters[1]]

    # gap penalty d
    d = parameters[2]

    # save all possible two sequences to align in an align object array
    all_alignments = []
    for i in range(len(sequences)):
        for j in range(len(sequences)-i-1):
            new_alignment = alignment(sequences[i], sequences[len(sequences)-j-1], s, d)
            all_alignments.append(new_alignment)

    # finds A*max and A*rest
    max_rest = find_max_and_rest(all_alignments)
    
    print(max_rest[1].first)
    print(max_rest[1].second)
    print(max_rest[0].first)
    print(max_rest[0].second)

    # finds A*cross
    cross = find_cross(max_rest)

    # T3




if __name__ == "__main__":
    # try:
        args = create_parser()
        # accesing the path of the files
        print(args.file_one)

        main()
    # except:
        # print('Try:  python3 template-A3.py -f1 to_msa.fasta')
