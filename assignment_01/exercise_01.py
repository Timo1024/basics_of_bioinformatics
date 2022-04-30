import argparse
import math
import contextlib

def check_if_valid_chars(char_array, string, linebreak = True, space = True):
    '''
    checks if the string is made of the alphabet char_array
    '''

    # if linebreaks at the end are allowed then add "\n" to the alphabet
    if(linebreak):
        char_array.append("\n")

    if(space):
        char_array.append(" ")

    # checks if every char in the string is correct
    for char in string:
        if not char in char_array:
            print(char)
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

def write_new_file(file, filename):
    '''
    uses the given file to get the headers
    then the user can set the sequences for these headers in the command line
    then the file is saved as a new file named filename
    '''

    # saves all headings
    heading_array = []
    with open(file) as file_content:
        for line in file_content:
            if(line.startswith(">")):
                heading_array.append(line)

    # creates new File
    f = open(filename, "w")

    # user can type the new sequences and the heading and the new sequences
    # are added to the new file
    for heading in heading_array:
        f.write(heading)
        cons_output = "Type the new Sequence for"
        cons_output += heading

        # checks if input is a correct sequence
        correct_input = False
        while not correct_input:
            new_sequence = input(cons_output)
            correct_input = check_if_valid_chars(["A", "T", "G", "C"], new_sequence)
        new_sequence += "\n"
        f.write(new_sequence)

    f.close()

    print("your new file was saved as", filename)

def make_reverse_complement(sequence):
    '''
    makes the reverse complement of a sequence
    '''

    # reverse the sequence
    sequence = sequence[::-1]

    # changes every letter with the reverse letter
    new_sequence = ""
    for base in sequence:
        if(base == "A"):
            new_sequence += "T"
        elif(base == "T"):
            new_sequence += "A"
        elif(base == "G"):
            new_sequence += "C"
        else:
            new_sequence += "G"

    return new_sequence

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


def reverse_sequence_order_and_save_reverse_comlements(file, filename):
    '''
    reverses the order of the sequences and saves the reverse complement of the
    sequences
    '''

    # saves all headings and sequences in two arrays
    heading_sequence_array = extract_headings_sequences(file)
    heading_array  = heading_sequence_array[0]
    sequence_array = heading_sequence_array[1]

    # makes the reverse complement of all sequences
    reverse_complement_array = []
    for sequence in sequence_array:
        sequence = "".join(sequence)
        reverse_complement_array.append(make_reverse_complement(sequence))

    # creates new File
    f = open(filename, "w")

    # reverse heading_array and reverse_complement_array
    heading_array = heading_array[::-1]
    reverse_complement_array = reverse_complement_array[::-1]

    # adding headings and sequences to new file
    count = 0
    for heading in heading_array:
        heading += "\n"
        f.write(heading)
        sequence = reverse_complement_array[count]
        sequence = '\n'.join(sequence[i:i+60] for i in range(0, len(sequence), 60))
        sequence += "\n"
        f.write(sequence)
        count += 1

    f.close()

    print("the reverse complement of the sequences was saved as", filename)

def calc_sp(sequence_array, aa_array):
    '''
    calculates the single propabilities of all amino acids
    and returns them in an array (positions are the same as in the aa_array)
    '''

    # concat all the sequences to one string
    sequence = ""
    for s in sequence_array:
        sequence += s[0]
    
    # save length of all sequences combined
    length = len(sequence)

    # initialize array with all sps
    sp_array = []

    for aa in aa_array:
        count = 0
        for char in sequence:
            if(char == aa):
                count += 1
        sp_array.append(count/length)
        
    return sp_array


def how_often_match_x_y(char_to_compare_array, char_1, char_2):
    '''
    calculates how often the two chars match in the given array
    '''

    # counts the matches
    match_count = 0

    for aa in char_to_compare_array:
        if aa == char_1:

            # subtract 1 b/c it counts itself as match when both chars are equal
            if(char_1 == char_2):
                match_count -= 1

            for bb in char_to_compare_array:
                if bb == char_2:
                    match_count += 1

    # if the two chars are equal every match is counted twice, 
    # so it has to be devided by 2
    if char_1 == char_2:
        match_count = match_count/2

    return match_count


def calc_jp(sequence_array, aa_array):
    '''
    calculates the joint propabilities of all amino acid pairs
    and returns them in an matrix (positions are the same as in the aa_array)
    '''

    # initialize matrix to store joint propabilities
    jp_matrix = []

    # calculate how many possible alignments are possible
    number_sequences = len(sequence_array)
    length_sequence = len(sequence_array[0][0])
    possible_alignments = sum(list(range(number_sequences))) * length_sequence
    
    # calculate jp for every combination
    # TODO it is enough to calculate one half of the matrix and then mirror it
    for aa1 in aa_array:
        jp_matrix.append([])
        for aa2 in aa_array:
            matches = 0
            for i in range(length_sequence):
                char_to_compare_array = []
                for array in sequence_array:
                    char_to_compare_array.append(array[0][i:i+1])
                matches += how_often_match_x_y(char_to_compare_array, aa1, aa2)
            jp = matches / possible_alignments

            # add jp to matrix
            jp_matrix[aa_array[aa1]].append(jp)

    return jp_matrix

def calc_substitution_matrix(single_propability_array, joint_propability_matrix, aa_array, no_match_score):
    '''
    calculates the substitution matrix from the given single propabilities ans joint propabilities
    '''

    # initialize substitution matrix
    sm = []

    for aa1 in aa_array:
        sm.append([])
        for aa2 in aa_array:
            if(((single_propability_array[aa_array[aa1]] * single_propability_array[aa_array[aa2]]) == 0) or
                (joint_propability_matrix[aa_array[aa1]][aa_array[aa2]] == 0)):
                sm[aa_array[aa1]].append(no_match_score)
            else:
                sm[aa_array[aa1]].append(
                    math.log(
                        joint_propability_matrix[aa_array[aa1]][aa_array[aa2]] / 
                        (single_propability_array[aa_array[aa1]] * single_propability_array[aa_array[aa2]]),
                        2
                    )
                )

    return sm


def print_matrix(substitution_matrix, aa_array, txt = False):
    '''
    prints the substitution Matrix in the command line
    '''

    # prints first two lines of matrix with all the amino acids
    if txt:
        print("\t /  \t", end = '')
        for char in aa_array:
            print("\t    ", end = '')
    else:
        print("\t /\t", end = '')
        for char in aa_array:
            print("\t", end = '')
    print("\\")

    print("\t|   \t", end = '')
    for char in aa_array:
        print(char, "  \t", end = '')
    print(" |")

    # calculate middle of the height of the matrix
    middle = math.floor((len(aa_array)-1)/2)

    counter = 0
    for char1 in aa_array:
        # print one line at a time
        if(counter == middle):
            print("S=\t| ", char1, "\t", end = '')
        else:
            print("\t| ", char1, "\t", end = '')
        for char2 in aa_array:
            value = round(substitution_matrix[aa_array[char1]][aa_array[char2]], 2)
            print(value, " \t", end = '')
        print(" |")
        counter += 1

    # print last line of matrix
    if txt:
        print("\t \\  \t", end = '')
        for char in aa_array:
            print("\t    ", end = '')
    else:
        print("\t \\\t", end = '')
        for char in aa_array:
            print("\t", end = '')
    print("/")



def compute_substitution_matrix(file, no_match_score, txt = False):
    '''
    computes the substitution matrix for the 
    sequences given in the file handed to the function
    prints out the matrix in the command line
    '''

    # get the sequences in one array
    heading_sequence_array = extract_headings_sequences(file)
    sequence_array = heading_sequence_array[1]

    # makes an associative array
    # the number represents the index, where the amino acid is located in the arrays
    aa_array = {
        "A": 0, "R": 1, "N": 2, "D": 3, 
        "C": 4, "E": 5, "Q": 6, "G": 7, 
        "H": 8, "I": 9, "L": 10, "K": 11, 
        "M": 12, "F": 13, "P": 14, "S": 15, 
        "T": 16, "W": 17, "Y": 18, "V": 19
    }
    # TODO delete commented stuff below
    # aa_array = {
    #     "A": 0, "R": 1, "C": 2
    # }

    # calculate the single propability for every aa
    single_propability_array = calc_sp(sequence_array, aa_array)

    # calculate the joined propability for every aa pair
    joint_propability_matrix = calc_jp(sequence_array, aa_array)

    # calculate substitution matrix
    substitution_matrix = calc_substitution_matrix(single_propability_array, joint_propability_matrix, aa_array, no_match_score)

    # print substitution matrix
    print_matrix(substitution_matrix, aa_array, txt)


def print_heading(string):
    '''
    prints a styled heading with the given string
    '''

    print("")
    print("--------------------------------------------------")
    print(string)
    print("--------------------------------------------------")


def main():
    '''
    The main function should contain all functions that solved the stated tasks. 
    '''
    # T2.a
    # TODO schauen, was wirklich in die main Funktion rein soll (vermutlich nur der Teil bei a))

    new_file_name = "MultipleSeqsReverse.fasta"
    file1 = args.file_one

    # all base chars
    base_chars = ["A", "T", "G", "C"]

    if(check_if_fasta(file1, base_chars)):
        
        # counting the bases of the sequences
        print_heading("Counting of the nucleotides in the sequences:")
        count_nucleotides(file1)

        # calculating the reverse complement of the sequences
        print_heading("Calculating the reverse complement of the sequences:")
        if(not new_file_name.endswith(".fasta")):
            new_file_name += ".fasta"
        reverse_sequence_order_and_save_reverse_comlements(file1, new_file_name)

    else:
        print("wrong file format. Please insert a .fasta file with bases")

    # T3

    file2 = args.file_two

    # all amino acid chars
    amino_chars = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    if(check_if_fasta(file2, amino_chars)):

        # calculating the substitution matrix
        # it needs a score for the case when the jp or sp is equal to zero
        print_heading("Substitution matrix:")
        compute_substitution_matrix(file2, -10)

        # save substitution matrix in txt file
        file_path = 'substitution_matrix.txt'
        with open(file_path, "w") as o:
            with contextlib.redirect_stdout(o):
                compute_substitution_matrix(file2, -10, True)


    else:
        print("wrong file format. Please insert a .fasta file with amino acids")




if __name__ == "__main__":
    try:
        args = create_parser()
        # accesing the path of the files
        print(args.file_one)
        print(args.file_two)

        main()
    except:
        print('Try:  python3 exercise_01.py -f1 MultipleSeqs.fasta -f2 msa-scoring-matrix.fasta')
