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
        matrix[i][i-1] = [0,[[0,0]]]

    for i in range(L):
        matrix[i][i] = [0,[[0,0]]]

    for n in range(1, L):
        for j in range(n, L):

            i = j-n

            if(j - i > l):

                # calculate maximum
                if(a[i] + a[j] in base_pairs):
                    delta = 1
                else: 
                    delta = 0

                fourth_entry_array = [0]
                for k in range(i+1, j):
                    fourth_entry_array.append(matrix[i][k][0] + matrix[k+1][j][0])

                max_fourth_entry = max(fourth_entry_array)
                indices_fourth_entry = list(np.where(np.array(fourth_entry_array) == max_fourth_entry)[0])

                traceback_fourth_entry = []
                for index in indices_fourth_entry:
                    traceback_fourth_entry.append([[i, i+index], [i+index+1, j]])

                maximum_array = [
                    matrix[i+1][j][0], 
                    matrix[i][j-1][0], 
                    matrix[i+1][j-1][0] + delta,
                    max_fourth_entry
                ]

                maximum = max(maximum_array)
                matrix[i][j][0] = maximum

                # fill Traceback in the matrix
                indices = list(np.where(np.array(maximum_array) == maximum)[0])
                traceback = []
                for index in indices:
                    if(index == 0):
                        traceback.append([[i+1, j]])
                    elif(index == 1):
                        traceback.append([[i, j-1]])
                    elif(index == 2):
                        traceback.append([[i+1, j-1]])
                    else: # index == 3
                        traceback += traceback_fourth_entry
                matrix[i][j][1] = traceback[0] # TODO delete [0] for all tracebacks
            else:
                matrix[i][j][0] = 0
                matrix[i][j][1] = [[0, 0]]
             
    return matrix


def traceback(matrix):
    """
    returns all nucleotide pairs for one traceback
    """
    n = len(matrix)
    tracebacks = matrix[0][n-1][1]
    pairs = []

    while(tracebacks[0] != [0, 0]):
        new_tracebacks = []
        for i in tracebacks:
            
            if(
                len(matrix[i[0]][i[1]][1]) == 1 and 
                matrix[i[0]][i[1]][1][0][0] == i[0]+1 and
                matrix[i[0]][i[1]][1][0][1] == i[1]-1
            ):
                pairs.append([i[0], i[1]])
                
            new_tracebacks += matrix[i[0]][i[1]][1]
        tracebacks = new_tracebacks

    return pairs


def make_dot_bracket(pairs, length):
    """
    input: nucleotide pairs
    output: dot-bracket notation
    """
    dot_bracket = ["."]*length
    for i in pairs:
        dot_bracket[min(i)] = "("
        dot_bracket[max(i)] = ")"
    return "".join(dot_bracket) 


def main():
    '''
    The main function
    '''

    file        = args.file
    loop_length = args.length

    headings_sequences = extract_headings_sequences(file)

    for heading_sequence_index in range(len(headings_sequences[0])):
        
        print("")

        heading  = headings_sequences[0][heading_sequence_index]
        sequence = headings_sequences[1][heading_sequence_index]
        print(heading)
        print(sequence)

        dp_matrix = compute(sequence, loop_length)

        length = len(dp_matrix)
        pairs = traceback(dp_matrix)
        dot_brackets = make_dot_bracket(pairs, length)

        print(dot_brackets)

        # get maximum amount of base pairs
        amount_basepairs = dp_matrix[0][length-1][0]
        print("The maximum number of basepairs is " + str(amount_basepairs))

    print("")


if __name__ == "__main__":
    try:
        args = create_parser()
        
        # accessing the arguments of argparse
        print("fasta file: " + args.file)
        print("minimal length for a loop: " + str(args.length))

        main()
    except:
        print('Try:  python3 Sara_Kemmler_Robin_Bonkass_A11.py -f nussinov.fasta -l 3')