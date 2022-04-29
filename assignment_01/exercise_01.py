import argparse
from asyncio.windows_events import NULL
from itertools import count

def count_nucleotides(file):
    '''
    prints out the length of all nucleotide sequences in the given fasta file
    '''

    with open(file) as file_content:

        # initialize varibles
        count_nucleotides = 0
        name_of_sequence = ""

        # loops through all lines in the given file
        for line in file_content:

            # a ">" at the start of a line indicates that there is a new sequence
            # every next start of the next sequence the program prints the result of the last sequence
            if(line.startswith('>')):
                if((count_nucleotides != 0) and (name_of_sequence != "")):
                    print("length of", name_of_sequence[:name_of_sequence.rfind("\n")], "=", count_nucleotides)
                    
                    # resetting the count of nucleotides for the next sequence
                    count_nucleotides = 0
                
                # saves the name of the next sequence for later, when ist printed out
                name_of_sequence = line

            # if the line doesnt start with a ">", the program counts the nucleotides in the line
            # add adds them to the nucleotide count   
            if(line.startswith('>') == False):
                for k in line:

                    # checks if the current char is one of the nucleotides
                    # this prevents the linebreak to count as a nucleotide
                    if k in ["A", "T", "G", "C"]:
                        count_nucleotides = count_nucleotides + 1

        # prints the length of the last nucleotide sequence
        if((count_nucleotides != 0) and (name_of_sequence != "")):
            print("length of", name_of_sequence[:name_of_sequence.rfind("\n")], "=", count_nucleotides)
            count_nucleotides = 0

def create_parser():
    '''
    Creates a argument parser to facilitate file reading.

    '''
    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-f1', '--file-one',
                   help="File for the the first programming task")
    p.add_argument('-f2', '--file-two',
                   help="File for the the second programming task")

    return(p.parse_args())


def main():
    '''
    The main function should contain all functions that solved the stated tasks. 
    '''
    # T2.a

    count_nucleotides(args.file_one)

    # TODO T3




if __name__ == "__main__":
    # try:
        args = create_parser()
        # accesing the path of the files
        print(args.file_one)
        print(args.file_two)

        main()
    # except:
    #     print('Try:  python3 template-a1.py -f1 MultipleSeqs.fasta -f2 msa-scoring-matrix.fasta')
