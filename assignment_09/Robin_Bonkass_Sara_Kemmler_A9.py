from ctypes import Array
from operator import indexOf
import numpy as np
import argparse

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
    return -np.inf if n == 0 else np.log(n)

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
        # print("Number of states: \t" + str(self.state_number))

        # Name of states
        self.state_names = clean_array(lines[3].split(" "))
        # print("Names of states: \t" + str(self.state_names))

        # Number of Symbols
        self.symbols_number = int(lines[5])
        # print("Number of symbols: \t" + str(self.symbols_number))

        # Name of symbols
        self.symbol_names = clean_array(lines[7].split(" "))
        # print("Names of symbols: \t" + str(self.symbol_names))

        # Transition probabilities matrix
        list_matrix_transition = []
        # for i in range(9, 18):
        for i in range(9, 12):
            list_matrix_transition.append(list(map(float, clean_array(lines[i].split(" ")))))
        self.transm_matrix = np.matrix(list_matrix_transition)
        # print("Transition Matrix:")
        # print(self.transm_matrix)

        # Emission probabilities matrix
        list_matrix_emission = []
        # for i in range(19, 28):
        for i in range(13, 16):
            list_matrix_emission.append(list(map(float, clean_array(lines[i].split(" ")))))
        self.emission_matrix = np.matrix(list_matrix_emission)
        # print("Emission Matrix:")
        # print(self.emission_matrix)
            


    def get_transition_probability(self, start_state, to_state):
        """
        Optional function: Given states from the current HMM, returns the transition probability.

        Args:
            start_state (str): starting state
            to_state (str): ending state

        Returns:
            float: probability of going from start_state to to_state
        """
        pass

    def get_emission_probability(self, state, symbol):
        """
        Optional function: Given a state and a symbol from the current HMM, returns the emission probability of the symbol under the state.

        Args:
            state (str): current state
            symbol (str): symbol that is emitted

        Returns:
            float: probability of emitting symbol under state
        """
        pass

    def runViterbi(self, sequence):
        """ For the given sequence of symbols, return the decoded sequence of symbols using the Viterbi algorithm.

        Args:
            sequence (str): a string of symbols to decode using the Viterbi algorithm.

        Returns:
            str: decoded sequence of states for the given sequence.
        """

        length = len(sequence)

        decoded_states = ""
        # TODO: implement the following functions.
        # You may need to adapt the input parameters for each function.

        # viterbi_init: initialization of the Viterbi algorithm
        dp_matrix = self.viterbi_init(length)

        # viterbi_matrix_filling: fills up the matrices
        dp_matrix = self.viterbi_matrix_filling(dp_matrix, sequence)

        # viterbi_terminate: computes the final step of the recursion.
        # viterbi_traceback: returns the decoded sequence of states for the given sequence of symbols

        return decoded_states

    def viterbi_init(self, length):

        dp_matrix = np.empty((self.state_number, length+2))
        dp_matrix[0, 0] = 1
        for i in range(1, self.state_number):
            dp_matrix[i, 0] = 0

        return dp_matrix

    def viterbi_matrix_filling(self, dp_matrix, sequence):

        # e_l(x_i) * max_k(v_k(i-1) * p_kl)

        # xi are the letters in the sequence
        # l are the states

        for i in range(len(sequence)): # i is index in the max part of the expression
            for l in range(self.state_number): # l is the index of the state

                # calculate the e_l(x_i)
                elxi = self.emission_matrix[l, indexOf(self.symbol_names, sequence[i])]
                
                max_k_array = []
                for k in range(self.state_number):
                    print(l)
                    max_k_array.append(dp_matrix[k, i] * self.transm_matrix[k, l-1])
                # print(max(max_k_array))
                dp_matrix[l, i+1] = elxi * max(max_k_array)

                # print(str(elxi) + " * max(" + str(max_k_array) + ")")

                # print(dp_matrix)

                pass

    def viterbi_terminate(self):
        pass

    def viterbi_get_path(self):
        pass


def prettyPrinting(input, decoded):
    """ Formats the given sequences to match the desired output format.

    Args:
        input (str): the original input sequences
        decoded (str): the sequence with the decoded states
    """
    pass


def main():
    hmm_object = HMMhandler()
    hmm_object.read_hmm(args.hmm_model)

    # TODO Parse fasta file for sequences
    original = extract_headings_sequences(args.sequences)[1]

    # TODO For each sequence in the fasta file run the viterbi algorithm.
    decoded = []
    for sequence in original:
        decoded.append(hmm_object.runViterbi(sequence))
    print(decoded)
    
    # TODO Once decoded, print the original and the decoded sequences with the desired output format.


if __name__ == "__main__":
    # try:
        args = create_parser()
        main()
    # except:
    #     print('Try:  python3 Robin_Bonkass_Sara_Kemmler_A9.py -seqs input_hmm.fasta -hmm cpg.hmm')
