"""

Bioinformatics
Assignment 3
Excercise 4: Feng Doolittle Distance

Robin Bonkaß
Sara Kemmler

TODO description of code

"""

# TODO delete in all files unused imports

from ctypes import alignment
from fd_lib import *

def main():

    """
    Task 3
    """

    '''
    initialize all needed values
    '''

    # extracting the sequences of the file into one array
    file1 = args.file_one
    sequences = extract_headings_sequences(file1)[1]

    # reads parameters for the alignment
    parameters = read_arguments(args)

    # scoring matrix [match score, mismatch score]
    s = [parameters[0], parameters[1]]

    # gap penalty d
    d = parameters[2]

    matrix = calc_distance_matrix(sequences, s, d)

    # prints results in distance_matrix.txt
    make_file(matrix, s, d, file1 , sequences)

if __name__ == "__main__":
    try:
        args = create_parser()

        # accesing the path of the files
        print(args.file_one)

        main()
    except:
        print('Try:  python3 feng_doolittle.py -sAA 3 -sAB -2 -d 4 -f1 to_msa.fasta')