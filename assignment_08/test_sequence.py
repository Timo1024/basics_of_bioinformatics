'''
Assignment 8
Task 4b)
Robin Bonkaß
Sara Kemmler

This program gets a fasta file with a sequence
and a P+ and P- Matrix as files
The program prints out into the console the 
probabilities of the sequence under the two 
given models and the respective log-odds ratio
'''

import math


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


def get_matrix_from_file(transition_matrix_file):
    '''
    returns matrix as array from given file
    '''

    matrix = []

    with open(transition_matrix_file) as file_content:
        for line in file_content:
            if(line[0:2] == "0."):

                new_line = line.split(" ")
                if(new_line[-1:][0] == "\n"):
                    new_line = new_line[:-1]

                for i in range(len(new_line)):
                    new_line[i] = float(new_line[i])
                
                matrix.append(new_line)
    
    return matrix


def compute_probability(transition_matrix_file, contig_file):
    '''
    calculates probability for one contig for one transition matrix
    '''

    sequence = "b" + extract_headings_sequences(contig_file)[1][0] + "e"

    matrix = get_matrix_from_file(transition_matrix_file)

    transitions = {
        'bb' : [4,4],
        'bA' : [4,0],
        'bG' : [4,2],
        'bC' : [4,1],
        'bT' : [4,3],
        'be' : [4,5],
        'Ab' : [4,4],
        'AA' : [0,0],
        'AG' : [0,2],
        'AC' : [0,1],
        'AT' : [0,3],
        'Ae' : [0,5],
        'Gb' : [4,4],
        'GA' : [2,0],
        'GG' : [2,2],
        'GC' : [2,1],
        'GT' : [2,3],
        'Ge' : [2,5],
        'Cb' : [4,4],
        'CA' : [1,0],
        'CG' : [1,2],
        'CC' : [1,1],
        'CT' : [1,3],
        'Ce' : [1,5],
        'Tb' : [4,4],
        'TA' : [3,0],
        'TG' : [3,2],
        'TC' : [3,1],
        'TT' : [3,3],
        'Te' : [3,5],
        'eb' : [4,4],
        'eA' : [5,0],
        'eG' : [5,2],
        'eC' : [5,1],
        'eT' : [5,3],
        'ee' : [5,5],
    }

    probability = -1

    for i in range(len(sequence)-1):
        if(matrix[transitions[sequence[i+0:i+2]][0]][transitions[sequence[i+0:i+2]][1]] != 0):
            if(probability == -1):
                probability = matrix[transitions[sequence[i+0:i+2]][0]][transitions[sequence[i+0:i+2]][1]]
            else:
                probability *= matrix[transitions[sequence[i+0:i+2]][0]][transitions[sequence[i+0:i+2]][1]]

    return probability


def compute_log_odds(transition_matrix_file_plus, transition_matrix_file_minus, contig_file):
    '''
    computes the log odds for a contig
    '''

    return math.log(compute_probability(transition_matrix_file_plus, contig_file)/compute_probability(transition_matrix_file_minus, contig_file))