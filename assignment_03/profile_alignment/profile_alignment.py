"""

Bioinformatics
Assignment 3
Excercise 3: Profile Alignment

Robin Bonkaß
Sara Kemmler

you can give a .fasta file with multiple sequences via command in the command 
line. The program extracts all sequences from the given file. It allways
searches for the two best fitting sequences and aligns them or their existing
profile alignments. Then it uses the next two best fitting ones. After all
of the sequences were aligned one with each other the program starts again. But
now with the calculated profile alignments. The two sequences used to align two
profile alignments are the sequences which get the best score when aligned
pairwise.

The program is functional for as many sequences as you provide. Especially, the
number of alignments doesn't have to be a power of 2. Because if there is one
sequence (at the start) or one profile alignment left over, which can't be
aligned with another one, it just gets used in the next iteration as it is.

The program finishes if just one profile alignment is left over. It returns
the resulting alignment in the command line and prints it to a .txt file
named profile_alignment.txt.

To start the program type:
python3 profile_alignment.py -sAA 3 -sAB -2 -d 4 -f1 to_msa.fasta

"""


from pa_lib import *

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

    paList = initialize_paList(sequences)

    '''
    loop through the building of the new pa list until there is just one pa left
    '''

    while (len(paList) > 1):
        new_paList = []

        '''
        loop through the current pa list and calculate the highest s_cross
        and then aligning the two pa's with the highest s_corss
        until the pa list is empty and then use the new pa list as the pa list
        '''

        while(len(paList) > 1):

            result = calc_scross(paList, s, d)

            # update paList b/c two pas were deleted
            paList = result[0]

            align = pair_alignment(result[1][1], result[2][1], s, d)

            pa1 = result[1][0]
            pa2 = result[2][0]

            new_pa_array = align.getAlignmentPa(pa1, pa2)[2]

            new_pa = profile_alignment(new_pa_array)

            new_paList.append(new_pa)

        # if there was an uneven number of pa's in the list, the last one
        # gets no alignment and is just added to the new pa list as it is
        if(len(paList) == 1):
            new_paList.append(paList[0])

        paList = new_paList

    # prints the sequences of the last pa in the command line
    # so that you can see how the profile alignment looks like
    print("")
    for i in paList[0].seqs:
        print(i)
    print("")

    # prints results in profile_alignment.txt
    make_file(paList[0].seqs, s, d, file1, sequences)


if __name__ == "__main__":
    try:
        args = create_parser()

        # accesing the path of the files
        print(args.file_one)

        main()
    except:
        print('Try:  python3 profile_alignment.py -sAA 3 -sAB -2 -d 4 -f1 to_msa.fasta')

