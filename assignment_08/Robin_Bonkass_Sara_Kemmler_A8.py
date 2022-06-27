'''
Assignment 8
Robin Bonkaß
Sara Kemmler

runs the program
'''

from train_hmm import *
from test_sequence import *
import argparse

def create_parser():
    '''
    Creates a argument parser to facilitate file reading.
    '''

    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-cds', '--cds-file',
                   help="Training data for cds")

    p.add_argument('-notcds', '--not-cds-file',
                   help="Training data for not cds")

    p.add_argument('-contig', '--contig-file',
                   help="file containing one contig")

    p.add_argument('-namecds', '--name-cds-matrix-file',
                   help="the name for the file with the cds matrix")

    p.add_argument('-namenotcds', '--name-not-cds-matrix-file',
                   help="the name for the file with the not cds matrix")

    return(p.parse_args())


def main():
    '''
    runs the program
    '''

    cds = args.cds_file
    not_cds = args.not_cds_file
    contig = args.contig_file
    name_cds = args.name_cds_matrix_file
    name_not_cds = args.name_not_cds_matrix_file

    compute_transition_matrix(cds, name_cds)
    compute_transition_matrix(not_cds, name_not_cds)

    print("p(contig|cds) = ", end="")
    print(compute_probability(name_cds + ".txt", contig))
    print("p(contig|notcds) = ", end="")
    print(compute_probability(name_not_cds + ".txt", contig))
    print("log-odds = ", end="")
    print(compute_log_odds(name_cds + ".txt", name_not_cds + ".txt", contig))


if __name__ == "__main__":
    try:
        args = create_parser()

        main()
    except:
        print('Try:  python3 Robin_Bonkass_Sara_Kemmler_A8.py -cds ./material/cds_set.fasta -notcds ./material/notcds_set.fasta -contig ./material/contig.fasta -namecds cds_matrix -namenotcds notcds_matrix')

