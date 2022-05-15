from pa_lib import *

def main():
    """
    Task 2
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

    # TODO make nice looking output (maybe in new .txt file)



if __name__ == "__main__":
    try:
        args = create_parser()

        # accesing the path of the files
        print(args.file_one)

        main()
    except:
        print('Try:  python3 profile_alignment.py -f1 to_msa.fasta')

