import argparse
from fd_classes import *
import math

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

    p.add_argument('-sAA', '--score-AA',
                    help="score for matching pair")

    p.add_argument('-sAB', '--score-AB',
                    help="score for mismatching pair")

    p.add_argument('-d', '--gap-penalty',
                    help="score for the gap penalty")

    p.add_argument('-display', '--display',
                    help="0 = half matrix; 1 = half and diagonal entries; 2 = full matrix")

    return(p.parse_args())


def read_arguments(args):
    '''
    reads the match/mismatch and gap penalty from the command line
    '''
    return [int(args.score_AA), int(args.score_AB), int(args.gap_penalty), int(args.display)]


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

    print("----------------------")
    print("The distance Matrix D:")
    print("----------------------")
    
    for i in matrix:
        for j in i:
            if(j == "-"):
                print("-", end="")
            else:
                print(round(j, 4), end = "")
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
                f.write(str(round(j, 4)))
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


def calc_distance_matrix(sequences, s, d, display = 0):
    '''
    calculates the feng doolittle distance for every possible combination of
    the given sequences and saves them in an matrix
    '''
    matrix = []

    for i in range(len(sequences)):
        matrix.append([])
        for j in range(len(sequences)):
            matrix[i].append("-")

    substract = 1
    if(display == 1):
        substract = 0

    for i in range(len(sequences)):
        if(display != 2):
            for j in range(len(sequences)-i-substract):
                i_index = i
                j_index = len(sequences)-j-1
                new_alignment = pair_alignment(sequences[i_index], sequences[j_index], s, d)

                # calculating the distance score of the two sequences
                # Logarithm is with base e
                fd = -math.log(
                    (new_alignment.getSobs() - new_alignment.getSrand())/
                    (new_alignment.getSid() - new_alignment.getSrand()))

                matrix[i_index][j_index] = fd
        else:
            for j in range(len(sequences)):
                i_index = i
                j_index = len(sequences)-j-1
                new_alignment = pair_alignment(sequences[i_index], sequences[j_index], s, d)

                # calculating the distance score of the two sequences
                # Logarithm is with base e
                fd = -math.log(
                    (new_alignment.getSobs() - new_alignment.getSrand())/
                    (new_alignment.getSid() - new_alignment.getSrand()))

                matrix[i_index][j_index] = fd

    return matrix