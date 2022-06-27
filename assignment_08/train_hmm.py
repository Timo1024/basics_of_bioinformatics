'''
Assignment 8
Task 4a)
Robin Bonkaß
Sara Kemmler

This Program gets a .fasta file 
and a name for a new file and 
computes a transition matrix and 
saves it as <filename>.txt
'''


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


def init_transition_matrix():
    '''
    initializes the 6x6 transition matrix
    '''

    matrix = []
    for i in range(6):
        matrix.append([])
        for j in range(6):
            matrix[i].append(0)

    return matrix


def print_to_file(matrix, file_name):

    # makes file with the distance matrix
    f = open(file_name + '.txt', "w")

    f.write("# Generated transition matrix\n")
    f.write("# Number of states:\n")
    f.write("6\n")
    f.write("State labels: *=b, +=e\n")
    f.write("A C G T * +\n")
    f.write("# Transition matrix:\n")

    for i in matrix:
        for j in i:
            f.write(str(format(round(j, 4), ".4f")))
            f.write(" ")
        f.write("\n")  

    f.close()


def compute_transition_matrix(fasta_given, file_name):
    '''
    calculates and saves transition matrix
    '''

    sequences = extract_headings_sequences(fasta_given)[1]

    # calculating the observed frequencies transition matrix
    observed_transition_matrix = init_transition_matrix()

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
    
    for sequence in sequences:
        for i in range(len(sequence)):
            # for start
            if(i == 0):
                transition = transitions['b'+sequence[i]]
                observed_transition_matrix[transition[0]][transition[1]] += 1
            # for end
            if(i == len(sequence)-1):
                transition = transitions[sequence[i]+'e']
                observed_transition_matrix[transition[0]][transition[1]] += 1
            # for all x_i x_i-1
            if(i != 0):
                transition = transitions[sequence[i-1]+sequence[i]]
                observed_transition_matrix[transition[0]][transition[1]] += 1
        transition = transitions['ee']
        observed_transition_matrix[transition[0]][transition[1]] += 1

    # calculating the probabilities transition matrix
    probability_transition_matrix = init_transition_matrix()

    for i in range(len(probability_transition_matrix)):
        sum = 0
        for number in observed_transition_matrix[i]:
            sum += number
        for j in range(len(probability_transition_matrix[i])):
            try: 
                probability_transition_matrix[i][j] = observed_transition_matrix[i][j]/sum
            except:
                pass

    print_to_file(probability_transition_matrix, file_name)