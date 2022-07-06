
from math import exp
from operator import indexOf
import numpy as np
import sys


np.seterr(divide='ignore')

def read_fasta_file(file):
    '''
    reads fasta file and returns all sequences as a list
    '''
    sequences = []
    i = -1
    with open(file, "r") as f:
        for line in f:
            if line[:1] == ">":
                i += 1
                sequences.append("b")
            else:
                sequences[i] = sequences[i] + line.rstrip("\n")

    for i in range(0, len(sequences)):
        sequences[i] = sequences[i] + "e"

    #print(sequences)

    return(sequences)  


def log_data(n):
    """Optional function: Handles the case that a log recieves a 0. It returns minus infinity

    Args:
        n (float)

    Returns:
        float or minus inf: log of n or minus infinity if 0 is given
    """
    return -np.inf if n == 0 else np.log(n)


class HMMhandler():
    """A class for handling Hidden Markov Models
    """

    def __init__(self, state_names=[], symbol_names=[], transm_matrix=[], emission_matrix=[]):
        """Constructor for the HMM handler class. 
        We recommend saving the parameters of the HMM model as properties of this clas. 
        It is not a must.
        """
        self.state_number = len(state_names)
        self.symbols_number = len(symbol_names)
        self.state_names = state_names
        self.symbol_names = symbol_names
        self.transm_matrix = transm_matrix
        self.emission_matrix = emission_matrix

    def read_hmm(self, hmm_handler):
        """Parses an HMM from the fixed format. See script P. 163.
        Relies on the file having the same order and structure from the shown structure.
        return: list[number of states, states, number of symbols, Symbols, transition probabilities, emission probabilities]
        """
        parsed_List = []
        with open(hmm_handler, "r") as file:
            for line in file:
                if line[:1] == "#":
                    parsed_List = parsed_List
                else:
                    parsed_List.append(line.rstrip("\n"))
                    
        self.state_number = int(parsed_List[0])
        self.state_names = parsed_List[1].split(" ")
        self.symbol_number = int(parsed_List[2])
        self.symbol_names = parsed_List[3].split(" ")
        
        for i in range(4, len(parsed_List)):
            if i < 4 + self.state_number:
                self.transm_matrix.append(parsed_List[i].split(" "))
            if i >= 4 + self.state_number and i <= 5 + 2 * self.state_number:
                self.emission_matrix.append(parsed_List[i].split(" "))
        
        pass

    def get_transition_probability(self, start_state, to_state):
        """Optinal function: Given states from the current HMM, returns the transition probability.

        Args:
            start_state (str): starting state
            to_state (str): ending state

        Returns:
            float: probability of going from start_state to to_state
        """
        pass

    def get_emission_probability(self, state, symbol):
        """Optional function: Given a state and a symbol from the current HMM, returns the emission probability of the symbol under the state.

        Args:
            state (str): current state
            symbol (str): symbol that is emitted

        Returns:
            float: probability of emitting symbol under state
        """
        pass
    
    def index_to_state(self, traceback):
            '''
                Input: List of index states
                Output: States
            '''
            states = []
        
            for state_index in traceback[1:-1]:
                states.append(self.state_names[state_index])
            
            #print(states)
            return("".join(states))
        
        
    def runViterbi(self, sequence):
        """ For the given sequence of symbols, return the decoded sequence of symbols using the Viterbi algorithm.

        Args:
            sequence (str): a string of symbols to decode using the Viterbi algorithm.

        Returns:
            str: decoded sequence of states for the given sequence.
        """
        
        # TODO: implement the following functions.
        # You may need to adapt the input parameters for each function.
        
        self.sequence_as_list = list(sequence)
        self.sequence_length = len(sequence)
        decoded_states = ""
        self.viterbi_init()                                                 # viterbi_init: initializaion of the Viterbi algorithm
        self.viterbi_matrix_filling()                                       # viterbi_matrix_filling: fills up the matrices
        self.viterbi_terminate()                                            # viterbi_terminate: computes the final step of the recursion.
        self.viterbi_get_path()
        decoded_states = self.index_to_state(self.viterbi_get_path())       #viterbi_traceback: returns the decoded sequence of states for the given sequence of symbols
        

        return decoded_states

    def viterbi_init(self):
        
        self.virterbi_matrix = []
        
        for row in range(0, self.state_number):
            self.virterbi_matrix.append([])
            for column in range(0, self.sequence_length):
                if row == 0 and column == 0:
                    self.virterbi_matrix[row].append(float(1))
                else:
                    self.virterbi_matrix[row].append(float(0))
        
        #print(self.virterbi_matrix)
        pass
    
    def which_base(self, base):
        
        if self.sequence_as_list[base] == "A":
            return(0)
        if self.sequence_as_list[base] == "C":
            return(1)
        if self.sequence_as_list[base] == "G":
            return(2)
        if self.sequence_as_list[base] == "T":
            return(3)
    
    def viterbi_matrix_filling(self):

        counter = 1

        for base in self.sequence_as_list[1:-1]:
            for state_l in range(0, self.state_number):
                emission_base = indexOf(self.symbol_names, base)
                emission_probability = float(self.emission_matrix[state_l][emission_base])
                maximum_finder = []
                
                for state_k in range(0, self.state_number):
                    viterbi_k = float(self.virterbi_matrix[state_k][counter-1])
                    transmission = float(self.transm_matrix[state_k][state_l])
                    maximum_finder.append(np.exp(log_data(viterbi_k) + log_data(transmission)))
                
                maximum = max(maximum_finder)
                self.virterbi_matrix[state_l][counter] = emission_probability  * maximum 
                
            counter += 1

        self.print_matrix(self.virterbi_matrix)

    # FIXME only for debugging
    def print_matrix(self, matrix):
        for i in matrix:
            for k in i:
                if(len(str(np.round(k, 8))) > 3):
                    print(str(np.round(k, 8)) + "\t", end="")
                else:
                    print(str(np.round(k, 8)) + "\t\t", end="")
            print("")

    def viterbi_terminate(self):
        
        base = self.sequence_length - 2
        maximum_finder = []
        for state_k in range(0, self.state_number):
                    viterbi_k = float(self.virterbi_matrix[state_k][base])
                    transmission = float(self.transm_matrix[state_k][0])
                    maximum_finder.append(log_data(viterbi_k)+ log_data(transmission))
        
        self.virterbi_matrix[0][self.sequence_length-1] = max(maximum_finder)
        #print(self.virterbi_matrix)
        pass

    def viterbi_get_path(self):
        
        traceback_state = []
        
        for base in range(0, self.sequence_length):
            maximum = -np.Inf
            for state in range(0, self.state_number):
                if maximum <= self.virterbi_matrix[state][base]:
                    maximum = self.virterbi_matrix[state][base]
                    current_state = state
                    
            traceback_state.append(current_state)                

        return(traceback_state)

  


    def prettyPrinting(self, input, decoded):
        """ Formats the given sequences to match the desired output format.

        Args:
            input (str): the original input sequences
            decoded (str): the sequence with the decoded states
        """
        with open("output.txt", "w") as file:
            sequence = input[1:self.sequence_length-1]
            length = len(sequence)
            for i in range(0, length+1):
                if i * 60 > length + 1:
                    break
                file.write("".join(input[i*60:min(60*(i+1), length)]) + "\n")
                file.write("".join(decoded[i*60:min(60*(i+1), length)]) + "\n\n")
                
        pass


def main():
    hmm_object = HMMhandler()
    hmm_object.read_hmm("cpg.hmm")
    decoded = hmm_object.runViterbi("bTACGCAe")
    hmm_object.prettyPrinting("bACGCe", decoded)
    # TODO Parse fasta file for sequences
    # TODO For each sequence in the fasta file run the viterbi algorithm.
    # TODO Once decoded, print the original and the decoded sequences with the desired output format.


if __name__ == "__main__":
    main()
