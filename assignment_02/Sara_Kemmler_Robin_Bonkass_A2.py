import argparse
from asyncio.windows_events import NULL
from email.mime import base

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

def create_parser():
    '''
    Creates a argument parser to facilitate file reading.
    '''

    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-f1', '--file-one',
                   help="First file with sequence")
    p.add_argument('-f2', '--file-two',
                   help="Second file with sequence")

    return(p.parse_args())


def check_if_valid_chars(char_array, string, linebreak = True, space = True):
    '''
    checks if the string is made of the alphabet char_array
    '''

    # if linebreaks at the end are allowed then add "\n" to the alphabet
    if(linebreak):
        char_array.append("\n")

    # if spaces (occur sometimes at the end of a line) are allowed then add " " to the alphabet
    if(space):
        char_array.append(" ")

    # checks if every char in the string is correct
    for char in string:
        if not char in char_array:
            
            # prints, if a char is in file but not allowed
            print("Character", char, "is not allowed in this file")
            return False
    
    return True


def check_if_fasta(file, valid_chars):
    '''
    checks if the file is a valid .fasta file
    '''

    # init variables
    counter = 0
    old_line = ""

    with open(file) as file_content:
        for line in file_content:

            # checks whether the first line is a header
            if counter == 0:
                if(not line.startswith(">")):
                    return False

            # if the last line was a header the following line can't be a header
            if(old_line.startswith(">")):
                if(line.startswith(">")):
                    return False

            # check if the chars in a line of the sequence are valid
            if(not line.startswith(">")):
                if(not check_if_valid_chars(valid_chars, line)):
                    return False
            
            old_line = line
            counter = counter + 1

        # checks that last line is not a header
        if(old_line.startswith(">")):
            return False
    
    return True


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

    return [heading_array, sequence_array]


def compute(seq1, seq2, s, d):
    '''
    gets the two sequences seq1 and seq2, 
    the scoring matrix s with match and mismatch scores and
    the gap penalty score d
    computes the dynamic programming matrix with the Needleman-Wunsch algorithm
    returns the dynamic programming matrix
    in this matrix the traceback is included
    '''

    # length of first sequence
    n = len(seq1)

    # length of second sequnce
    m = len(seq2)

    # initialize matrix array
    matrix = []
    for j in range(m+1):
        matrix.append([])
        for i in range(n+1):
            matrix[j].append([0,[0,0]])

    # initialize matrix values
    for i in range(n+1):
            matrix[0][i] = [-i * d, [i-1, 0]]
    for j in range(m+1):
        matrix[j][0] = [-j * d, [0, j-1]]

    # iterates through all fields of the matrix
    for i in range(1,n+1):
        for j in range(1,m+1):

            # calculate if match or mismatch
            if(seq1[i-1] == seq2[j-1]):
                score = s[0]
            else:
                score = s[1]

            # calculate maximum
            maximum = max(matrix[j-1][i-1][0] + score, matrix[j-1][i][0]-d, matrix[j][i-1][0]-d)

            # sets the max value in the matrix and the traceback
            if(maximum == matrix[j-1][i-1][0] + score):
                matrix[j][i] = [maximum, [i-1,j-1]]
            elif(maximum == matrix[j-1][i][0]-d):
                matrix[j][i] = [maximum, [i,j-1]]
            else:
                # maximum == matrix[j][i-1][0]-d could be condition here
                matrix[j][i] = [maximum, [i-1,j]]

    return matrix


def traceback(matrix, seq1, seq2):
    '''
    calculates one alignment with the matrix
    and prints the optimal score
    '''

    print("the optimal score for an global alignment is ", matrix[len(matrix)-1][len(matrix[0])-1][0])

    # print_matrix(matrix, True)

    start_i = len(matrix)-1
    start_j = len(matrix[0])-1

    # print(matrix[5][8])

    seq1_aligned = ""
    seq2_aligned = ""

    old_i = len(matrix)-1
    old_j = len(matrix[0])-1

    count = 0

    while(start_j != -1 and start_i != -1):

        print("compare (", old_i, ",",  old_j , ") with (", start_i, ",",  start_j , ")" )
        
        if(old_i == start_i and (old_i != len(matrix)-1 or start_i != len(matrix)-1)):
            seq2_aligned += "-"
            seq1_aligned += seq1[old_j-1]
            # print(seq2_aligned)
        elif(old_j == start_j and (old_j != len(matrix[0])-1 or start_j != len(matrix[0])-1)):
            seq1_aligned += "-"
            seq2_aligned += seq2[old_i-1]
            # print(seq1_aligned)
        elif(count != 0):
            seq1_aligned += seq1[old_j-1]
            # print(seq2)
            # print(start_i-1)
            seq2_aligned += seq2[old_i-1]

        count += 1
        old_j = start_j
        old_i = start_i

        # print("j =", start_j, "; i =", start_i)

        # print_matrix(matrix)

        start_i = matrix[old_i][old_j][1][1]
        start_j = matrix[old_i][old_j][1][0]

    # print("j =", start_j, "; i =", start_i)

    print_matrix(matrix, False)
    print_matrix(matrix, True)

    print(seq1_aligned)
    print(seq2_aligned)


def main():
    '''
    The main function should contain all functions that solved the stated tasks. 
    '''

    base_chars = ["A", "T", "G", "C"]
    file1 = args.file_one
    file2 = args.file_two

    if(check_if_fasta(file1, base_chars) and (check_if_fasta(file2, base_chars))):
        
        # get the parameters from command line input
        parameters = read_arguments()

        # get the two sequences
        seq1 = "".join(extract_headings_sequences(file1)[1][0])
        seq2 = "".join(extract_headings_sequences(file2)[1][0])

        # scoring matrix [match score, mismatch score]
        s = [parameters[0], parameters[1]]

        # gap penalty d
        d = parameters[2]

        # computes the matrix
        matrix = compute(seq1, seq2, s, d)

        # calculates one alignment and the optimal score
        traceback(matrix, seq1, seq2)

    else:
        print("wrong files format. The sequences have to contain just nucleotides")


if __name__ == "__main__":
    # try:
        args = create_parser()
        
        # accesing the path of the files
        print(args.file_one)
        print(args.file_two)

        main()
    # except:
    #     print('Try:  python3 Sara_Kemmler_Robin_Bonkass_A2.py -f1 yersenia_1.fasta -f2 yersenia_2.fasta')