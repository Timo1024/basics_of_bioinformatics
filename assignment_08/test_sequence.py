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
            if(line[0:6] == "0.0000"):

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

# TODO delete
    print(probability)

    return probability


def compute_log_odds(transition_matrix_file_plus, transition_matrix_file_minus, contig_file):
    '''
    computes the log odds for a contig
    '''

    sequence = "b" + extract_headings_sequences(contig_file)[1][0] + "e"

    matrix_plus = get_matrix_from_file(transition_matrix_file_plus)
    matrix_minus = get_matrix_from_file(transition_matrix_file_minus)

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

    log_odds = 0

    for i in range(len(sequence)-1):
        cords = transitions[sequence[i:i+2]]
        if(matrix_minus[cords[0]][cords[1]] != 0):
            if(matrix_plus[cords[0]][cords[1]]/matrix_minus[cords[0]][cords[1]] != 0):
                log_odds += math.log(matrix_plus[cords[0]][cords[1]]/matrix_minus[cords[0]][cords[1]])

# TODO delete
    print(log_odds)

    return log_odds





# TODO delete or copy following to main .py file when finished

import argparse

def create_parser():
    '''
    Creates a argument parser to facilitate file reading.
    '''

    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-f1', '--file-one',
                   help="cds file with sequence")

    p.add_argument('-f2', '--file-two',
                   help="not cds file with sequence")

    p.add_argument('-f3', '--file-three',
                   help="contig file with sequence")

    return(p.parse_args())

def main():

    file1 = args.file_one
    file2 = args.file_two
    file3 = args.file_three

    print("p(x|p+)")
    compute_probability(file1, file3)
    print("p(x|p-)")
    compute_probability(file2, file3)
    print("log-odds")
    compute_log_odds(file1, file2, file3)

if __name__ == "__main__":
    # try:
        args = create_parser()
        
        # accessing the path of the files
        print(args.file_one)
        print(args.file_two)
        print(args.file_three)

        main()
    # except:
    #     print('Try:  python3 test_sequence.py -f1 tm1.txt -f2 tm2.txt -f3 ./material/contig.fasta')

