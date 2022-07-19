'''
Assignment 2
Robin Bonkaß
Sara Kemmler
'''

import argparse
import numpy as np

def create_parser():
    '''
    Creates a argument parser to facilitate file reading.
    '''

    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-f', '--file',
                   help="Fasta file containing sequences")
    p.add_argument('-l', '--length', type=int,
                   help="minimal length for a loop")

    return(p.parse_args())


def print_matrix(matrix, full = False):
    '''
    prints a matrix in command line
    used to debug
    '''
    
    for i in range(len(matrix)):
        if(not full):
            for j in range(len(matrix[0])):
                print(matrix[i][j][0], "\t", end = '')
            print("")
        else:
            for j in range(len(matrix[0])):
                if(j == 0 and i == 0):
                    print("\t", end = '')
                else:
                    print(matrix[i][j][1], "\t", end = '')
            print("")


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
                if(len(sequence_array[len(heading_array)-1]) > 0):
                    sequence_array[len(heading_array)-1] = ["".join(sequence_array[len(heading_array)-1])]

    for i in range(len(sequence_array)):
        sequence_array[i] = sequence_array[i][0]

    return [heading_array, sequence_array]


def compute(a, l):
    '''
    input:
    - the one sequences a,
    - a score for the minimal count of a loop

    processing:
    - computes the dynamic programming matrix with the Nussinov algorithm

    returns:
    - returns the dynamic programming matrix
      in this matrix the traceback is included
    '''

    print(a)

    # TODO traceback

    #  canonical base pairs
    base_pairs = ["AU", "UA", "GC", "CG"]

    # length of sequence
    L = len(a)

    # initialize matrix array with None
    matrix = []
    for i in range(L):
        matrix.append([])
        for j in range(L):
            matrix[i].append([None,None])

    # initialization of nussinov
    for i in range(1, L):
        matrix[i][i-1] = [0,[0,0]]

    for i in range(L):
        matrix[i][i] = [0,[0,0]]

    print_matrix(matrix, full = False)

    for n in range(1, L):
        for j in range(n, L):
            i = j-n

            # calculate maximum
            # print("i = " + str(i) + "; j = " + str(j) + "; n = " + str(n))

            if(a[i] + a[j] in base_pairs):
                delta = 1
            else: 
                delta = 0

            # print(a[i] + a[j] + " -> delta = " + str(delta))

            fourth_entry_array = [0]
            for k in range(i+1, j):
                fourth_entry_array.append(matrix[i][k][0] + matrix[k+1][j][0])

            # print("max i<k<j = " + str(fourth_entry_array))

            maximum_array = [
                matrix[i+1][j][0], 
                matrix[i][j-1][0], 
                matrix[i+1][j-1][0] + delta,
                max(fourth_entry_array)
            ]

            maximum = max(maximum_array)
            matrix[i][j][0] = maximum

            # fill Traceback in the matrix
            indices = list(np.where(np.array(maximum_array) == maximum)[0])
            print(indices)


    return matrix


def main():
    '''
    The main function
    '''

    file        = args.file
    loop_length = args.length

    headings_sequences = extract_headings_sequences(file)

    # TODO make nussinov for all sequences

    dp_matrix = compute(headings_sequences[1][0], loop_length)

    print_matrix(dp_matrix, False)



if __name__ == "__main__":
    # try:
        args = create_parser()
        
        # accesing the arguments of argparse
        print("fasta file: " + args.file)
        print("minimal length for a loop: " + str(args.length))

        main()
    # except:
    #     print('Try:  python3 Sara_Kemmler_Robin_Bonkass_A11.py -f nussinov.fasta -l 3')