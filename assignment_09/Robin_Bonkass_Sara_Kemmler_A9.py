from math import exp, log
from operator import indexOf
import numpy as np
import argparse
import decimal

# TODO überall Kommentare schreiben


np.seterr(divide='ignore')

def create_parser():
    '''
    Creates a argument parser to facilitate file reading.
    '''

    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-seqs', '--sequences',
                   help="The fasta file with the sequences")

    p.add_argument('-hmm', '--hmm-model',
                   help="File containing the HMM model")

    return(p.parse_args())

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


def log_data(n):
    """
    Optional function: Handles the case that a log recieves a 0. It returns minus infinity

    Args:
        n (float)

    Returns:
        float or minus inf: log of n or minus infinity if 0 is given
    """
    # return -np.inf if n == 0 else np.log(n)
    return decimal.Decimal('-Infinity') if n == 0 else decimal.Decimal(n).ln()

def clean_array(array):
    """
    returns a list without the last element if it is empty or a linebreak
    """
    if(array[-1:][0] == '' or array[-1:][0] == '\n' or array[-1:][0] == ' '):
        return array[:-1]
    else:
        return array


class HMMhandler():
    """
    A class for handling Hidden Markov Models
    """

    def __init__(self, state_names=[], symbol_names=[], transm_matrix=np.matrix([], []), emission_matrix=np.matrix([], [])):
        """
        Constructor for the HMM handler class. 
        We recommend saving the parameters of the HMM model as properties of this clas. 
        It is not a must.
        """
        self.state_number = len(state_names)
        self.symbols_number = len(symbol_names)
        self.state_names = state_names
        self.symbol_names = symbol_names
        self.transm_matrix = transm_matrix
        self.emission_matrix = emission_matrix

    def read_hmm(self, hmm_file):
        """
        Parses an HMM from the fixed format. See script P. 163.
        Relies on the file having the same order and structure from the shown structure.
        """

        lines = []
        with open(hmm_file) as file_content:
            for line in file_content:
                lines.append(clean_array(line))

        # number of states
        self.state_number = int(lines[1])

        # Name of states
        self.state_names = clean_array(lines[3].split(" "))

        # Number of Symbols
        self.symbols_number = int(lines[5])

        # Name of symbols
        self.symbol_names = clean_array(lines[7].split(" "))

        # Transition probabilities matrix
        list_matrix_transition = []
        for i in range(9, 9+self.state_number):
            list_matrix_transition.append(list(map(float, clean_array(lines[i].split(" ")))))
        self.transm_matrix = np.matrix(list_matrix_transition)

        # Emission probabilities matrix
        list_matrix_emission = []
        for i in range(10+self.state_number, 10+(self.state_number*2)):
            list_matrix_emission.append(list(map(float, clean_array(lines[i].split(" ")))))
        self.emission_matrix = np.matrix(list_matrix_emission)
            


    def get_transition_probability(self, start_state, to_state):
        """
        Given states from the current HMM, returns the transition probability.
        Important that the state e (0) is located in the first column of the transition matrix

        Args:
            start_state (str): starting state
            to_state (str): ending state

        Returns:
            float: probability of going from start_state to to_state
        """
        return self.transm_matrix[indexOf(self.state_names, start_state), indexOf(self.state_names, to_state)]

    def get_emission_probability(self, state, symbol):
        """
        Given a state and a symbol from the current HMM, returns the emission probability of the symbol under the state.

        Args:
            state (str): current state
            symbol (str): symbol that is emitted

        Returns:
            float: probability of emitting symbol under state
        """
        return self.emission_matrix[indexOf(self.state_names, state), indexOf(self.symbol_names, symbol)]

    def runViterbi(self, sequence):
        """ For the given sequence of symbols, return the decoded sequence of symbols using the Viterbi algorithm.

        Args:
            sequence (str): a string of symbols to decode using the Viterbi algorithm.

        Returns:
            str: decoded sequence of states for the given sequence.
        """

        length = len(sequence)

        # viterbi_init: initialization of the Viterbi algorithm
        dp_matrix = self.viterbi_init(length)

        # viterbi_matrix_filling: fills up the matrices
        dp_matrix = self.viterbi_matrix_filling(dp_matrix, sequence)

        # viterbi_terminate: computes the final step of the recursion.
        dp_matrix = self.viterbi_terminate(dp_matrix, sequence)

        # viterbi_traceback: returns the decoded sequence of states for the given sequence of symbols
        decoded_states = self.viterbi_get_path(dp_matrix, sequence)

        return decoded_states

    def viterbi_init(self, length):
        """
        initializes the dp matrix and the first column
        all elements of the matrix are a list containing two elements
        - the first is the value
        - the second is the row of the last column where the value came from
          (for traceback)
        """

        dp_matrix = np.empty((self.state_number, length+2), dtype=object)
        dp_matrix[0, 0] = [1, 0]
        for i in range(1, self.state_number):
            dp_matrix[i, 0] = [0, 0]

        return dp_matrix

    def viterbi_matrix_filling(self, dp_matrix, sequence):

        # e_l(x_i) * max_k(v_k(i-1) * p_kl)

        index_in_sequence = 0
        for xi in sequence: # i is the current symbol name
            for l in self.state_names: # l is the current state name

                # get emission probability of the symbol under the state
                elxi = self.get_emission_probability(l, xi)
                
                # get all values of the last column of the dp matrix times the transition probability
                max_k_array = []
                for k in range(self.state_number):
                    # max_k_array.append(dp_matrix[k, index_in_sequence][0] * self.get_transition_probability(self.state_names[k], l))
                    # max_k_array.append(exp(log_data(dp_matrix[k, index_in_sequence][0]) + log_data(self.get_transition_probability(self.state_names[k], l))))
                    max_k_array.append(decimal.Decimal(log_data(dp_matrix[k, index_in_sequence][0]) + log_data(self.get_transition_probability(self.state_names[k], l))).exp())
                # dp_matrix[indexOf(self.state_names, l), index_in_sequence+1] = [elxi * max(max_k_array), np.argmax(max_k_array)]
                arg_max = np.argmax(max_k_array)
                # arg_max = np.argmax(max_k_array) if np.argmax(max_k_array) != 0 else 1
                # dp_matrix[indexOf(self.state_names, l), index_in_sequence+1] = [exp(log_data(elxi) + log_data(max(max_k_array))), arg_max]
                # dp_matrix[indexOf(self.state_names, l), index_in_sequence+1] = [exp(log_data(elxi) + log_data(max(max_k_array))), arg_max]
                dp_matrix[indexOf(self.state_names, l), index_in_sequence+1] = [decimal.Decimal(log_data(elxi) + log_data(max(max_k_array))).exp(), arg_max]

                if(max(max_k_array) == 0 and l != "0"):
                    print("Error: max is 0 can't choose maximum")
                    print("Argmax is: " + str(np.argmax(max_k_array)))
                    print(max_k_array)
                    for k in range(self.state_number):
                        print("log(dp_matrix[k, index_in_sequence][0]) = " + str(log_data(dp_matrix[k, index_in_sequence][0])))
                        print("log(self.get_transition_probability(self.state_names[k], l)) = " + str(log_data(self.get_transition_probability(self.state_names[k], l))))
                        print("exp(v + transitionp) = " + str(log_data(dp_matrix[k, index_in_sequence][0]) + log_data(self.get_transition_probability(self.state_names[k], l))))
                    print("current index in sequence: " + str(index_in_sequence))
                    print("---------------------------------------------------")
                    raise NameError('max is 0')

            index_in_sequence += 1

        return dp_matrix

    def viterbi_terminate(self, dp_matrix, sequence):

        max_k_array = []
        for k in range(self.state_number):
            # max_k_array.append(dp_matrix[k, len(sequence)][0] * self.get_transition_probability(self.state_names[k], self.state_names[0]))
            # max_k_array.append(exp(log_data(dp_matrix[k, len(sequence)][0]) + log_data(self.get_transition_probability(self.state_names[k], self.state_names[0]))))
            max_k_array.append(decimal.Decimal(log_data(dp_matrix[k, len(sequence)][0]) + log_data(self.get_transition_probability(self.state_names[k], self.state_names[0]))).exp())
        dp_matrix[0, len(sequence)+1] = [max(max_k_array), np.argmax(max_k_array)]
        
        for i in range(1, self.state_number):
            dp_matrix[i, len(sequence)+1] = [0, 0]

        # print(dp_matrix)

        return dp_matrix
        

    def viterbi_get_path(self, dp_matrix, sequence):
        
        traceback = 0
        decoded = ""
        for i in reversed(range(len(sequence)+2)):
            decoded += self.state_names[traceback]
            traceback = dp_matrix[traceback, i][1]

        return decoded[::-1]
        



def prettyPrinting(input, decoded, f):
    """ Formats the given sequences to match the desired output format.

    Args:
        input (str): the original input sequences
        decoded (str): the sequence with the decoded states
    """
    # make arrays of length 60 of every string
    input_cut = []
    input_cut.append(input[0:50])
    for i in range(50, len(input), 60):
        input_cut.append(input[i:i+60])
    decoded_cut = []
    decoded_cut.append(decoded[0:51])
    for i in range(51, len(decoded), 60):
        decoded_cut.append(decoded[i:i+60])

    f.write("Symbols:  ")
    f.write(input_cut[0])
    f.write("\n")
    f.write("Viterbi: ")
    f.write(decoded_cut[0])
    f.write("\n\n")

    for i in range(1, len(input_cut)):
        f.write(input_cut[i])
        f.write("\n")
        f.write(decoded_cut[i])
        f.write("\n\n")


def main():

    print(str(exp(-743.7469247408213 + -1.7883043669741367)))

    hmm_object = HMMhandler()
    hmm_object.read_hmm(args.hmm_model)

    # Parse fasta file for sequences
    headings_sequences = extract_headings_sequences(args.sequences)
    original = headings_sequences[1]
    headings = headings_sequences[0]

    # For each sequence in the fasta file run the viterbi algorithm.
    decoded = []
    for sequence in original:
        decoded.append(hmm_object.runViterbi(sequence))
    
    # Once decoded, print the original and the decoded sequences with the desired output format.
    f = open("output.txt", "w")
    f.write("Original and decoded sequences after applying the viterbi algorithm:\n\n\n")

    for i in range(len(original)):
        f.write(headings[i])
        f.write("\n\n")
        prettyPrinting(original[i], decoded[i], f)
        f.write("\n")

    f.close()


if __name__ == "__main__":
    # try:
        args = create_parser()
        main()
    # except:
    #     print('Try:  python3 Robin_Bonkass_Sara_Kemmler_A9.py -seqs input_hmm.fasta -hmm cpg.hmm')
