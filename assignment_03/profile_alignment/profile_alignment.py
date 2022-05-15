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
    print(paList)


    while (len(paList) > 1):
        new_paList = []
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

            # print(paList)
            # print(new_paList[0].seqs[0])
            # print(new_paList[0].seqs[1])

        if(len(paList) == 1):
            new_paList.append(paList[0])

        paList = new_paList

    for i in paList[0].seqs:
        print(i)



if __name__ == "__main__":
    # try:
        args = create_parser()
        # accesing the path of the files
        print(args.file_one)

        main()
    # except:
        # TODO change file name and uncomment try-except
        # print('Try:  python3 template-A3.py -f1 to_msa.fasta')

