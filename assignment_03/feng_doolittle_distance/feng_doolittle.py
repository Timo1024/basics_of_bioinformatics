"""

Bioinformatics
Assignment 3
Excercise 4: Feng Doolittle Distance

Robin Bonkaß
Sara Kemmler

TODO description of code

TODO alle Argumente mit erstem Aufruf über geben (mit -sAA -sAB -d usw.)

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
    parameters = read_arguments()

    # scoring matrix [match score, mismatch score]
    s = [parameters[0], parameters[1]]

    # gap penalty d
    d = parameters[2]

    matrix = calc_distance_matrix(sequences, s, d)

    # prints results in profile_alignment.txt
    make_file(matrix, s, d, file1 , sequences)

if __name__ == "__main__":
    # TODO uncomment
    # try:
        args = create_parser()

        # accesing the path of the files
        print(args.file_one)

        main()
    # except:
    #     print('Try:  python3 feng_doolittle.py -f1 to_msa.fasta')