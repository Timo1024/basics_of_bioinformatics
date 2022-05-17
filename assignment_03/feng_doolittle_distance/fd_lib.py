import argparse
from fd_classes import *
import math
# TODO dürfen wir math library benutzen?

'''
Here are all functions which are used in the main file
'''

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


def make_file(matrix, s, d, file1, sequences):
    '''
    prints the distance matrix to the command line
    and to a file named file1
    '''
    # TODO print matrix to file

    print("----------------------")
    print("The distance Matrix D:")
    print("----------------------")
    
    for i in matrix:
        for j in i:
            if(j == "-"):
                print("-", end="")
            else:
                print(round(j, 3), end = "")
            print("\t", end = "")
        print("\n")

    # makes file with the distance matrix
    f = open("distance_matrix.txt", "w")

    heading = "------------------------------------------------------------------\nDistance Matrix using Feng-Doolittle distances for MSA\nused were "
    heading += str(len(sequences))
    heading += " Sequences from file "
    heading += file1
    f.write(heading)
    f.write("\n------------------------------------------------------------------\n\n")

    top_matrix = " /\t"
    for k in matrix[0]:
        top_matrix += "\t\t"
    top_matrix += "\\\n"
    f.write(top_matrix)

    for i in matrix:
        f.write("|\t")
        for j in i:
            if(j == "-"):
                f.write("-\t")
            else:
                f.write(str(round(j, 3)))
            f.write("\t")
        f.write(" |\n")

    bottom_matrix = " \\\t"
    for k in matrix[0]:
        bottom_matrix += "\t\t"
    bottom_matrix += "/\n"
    f.write(bottom_matrix)

    # writes the used algorithm parameters
    f.write("\nused Algorithm parameters:")
    f.write("\nMatch score = \t\t")
    f.write(str(s[0]))
    f.write("\nMismatch score = \t")
    f.write(str(s[1]))
    f.write("\nGap Penalty = \t\t")
    f.write(str(d))

    f.close()


def calc_distance_matrix(sequences, s, d):
    '''
    calculates the feng doolittle distance for every possible combination of
    the given sequences and saves them in an matrix
    '''
    matrix = []

    # TODO reicht es aus, wenn man eine Hälfte und nicht die Diagonale in der
    # Matrix ausfüllt?

    for i in range(len(sequences)):
        matrix.append([])
        for j in range(len(sequences)):
            matrix[i].append("-")


    for i in range(len(sequences)):
        for j in range(len(sequences)-i-1):
            i_index = i
            j_index = len(sequences)-j-1
            new_alignment = pair_alignment(sequences[i_index], sequences[j_index], s, d)

            # calculating the distance score of the two sequences
            fd = -math.log(
                (new_alignment.getSobs() - new_alignment.getSrand())/
                (new_alignment.getSid() - new_alignment.getSrand()), 10)

            matrix[i_index][j_index] = fd

    return matrix