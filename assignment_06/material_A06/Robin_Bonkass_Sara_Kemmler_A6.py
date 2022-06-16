import argparse
import math
import string
import random
import matplotlib.pyplot as plt
import numpy as np
import sys

def create_parser():
    '''
    Creates an argument parser to get arguments from command line.
    '''

    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-repeats', '--repeats', type=int,
                   help="Amount of repeats of the simulation for each initial size")

    p.add_argument('-consoleOutput', '--consoleOutput', action='store_true', 
                   help="makes output in command line showing the generations")
    p.add_argument('-fileOutput', '--fileOutput', dest='consoleOutput', action='store_false', 
                   help="print output in file \"Robin_Bonkass_Sara_Kemmler_generations.txt\"")
    p.set_defaults(consoleOutput=False)

    p.add_argument('-sizes', '--sizes', type=int, nargs="*",
                   help="n numbers which represent the initial population sizes")

    return(p.parse_args())

def get_random_pair(current_generation_array):
    """
    returns array with a "-" at one position instead of a letter
    -> this means that two letters are merged two one
    """

    length = len(current_generation_array)

    # counts all left letters
    letters_count = 0
    for i in range(length):
        if(current_generation_array[i] != "-"):
            letters_count += 1

    # get index for two letters which will be merged
    first_letter_index = math.floor(random.random() * letters_count)
    second_letter_index = math.floor(random.random() * letters_count)

    if(first_letter_index == second_letter_index):
        second_letter_index += 1
        second_letter_index = second_letter_index % letters_count

    # get first letter
    first = 0
    for i in range(length):
        if(first == first_letter_index):
            if(current_generation_array[i] != "-"):
                first_letter = current_generation_array[i]
                break
        else:
            if(current_generation_array[i] != "-"):
                first += 1

    # get second letter
    second = 0
    for i in range(length):
        if(second == second_letter_index):
            if(current_generation_array[i] != "-"):
                second_letter = current_generation_array[i]
                break
        else:
            if(current_generation_array[i] != "-"):
                second += 1

    # get lexicographic lower letter -> save this as first_letter and other as second_letter
    if(first_letter > second_letter):
        first_letter_temp = first_letter
        first_letter = second_letter
        second_letter = first_letter_temp

    # get position of second letter and make it a "-"
    merged_letter_index = current_generation_array.index(second_letter)
    current_generation_array[merged_letter_index] = "-"

    return current_generation_array

class PopGenSimulator:

    def __init__(self, size):
        if(size <= 0 or size > 26):
            print("Parameter size has to be an int between 1 and 26")
            raise NameError("wrong size for initial population")
        else:
            self.size = size  # population size
        
        self.time = 0

    def simulate_previous_generation(self, current_generation, k, waiting_time_overlap):
        """
        given a string of letters, each letter representing one individual,
        randomly generate previous generation. For this use use the exponential distribution exp(k(k-1)/2) for the waiting times.
        Record them.
        Then randomly choose a pair which coalesces.
        In previous generation, use lexicographically first letter of all
        children as label of parent, or -, if parent has no child
        :param current_generation: String
        :return: previous_generation
        """

        # make array from string
        current_generation_array=[]
        current_generation_array[:0]=current_generation

        # waiting_time = 0
        waiting_time_return = 0
        waiting_time = waiting_time_overlap

        if(waiting_time_overlap <= 1):
            if(waiting_time_overlap > 0):
                current_generation_array = get_random_pair(current_generation_array)
                k -= 1
            
            while(waiting_time<=1 and k >= 2):
                random_uniform = random.uniform(0,1)
                waiting_time += ((0-(math.comb(k,2) ** (-1))) * math.log(random_uniform)) * self.size
                self.time += ((0-(math.comb(k,2) ** (-1))) * math.log(random_uniform)) * self.size
                if(waiting_time <= 1):
                    current_generation_array = get_random_pair(current_generation_array)
                    k -= 1

        previous_generation = "".join(current_generation_array)
        
        if(waiting_time > 1):
            waiting_time_return = waiting_time - 1

        return [previous_generation, k, waiting_time_return]


    def is_MRCA_of_all_found(self, current_generation):
        """
        have we found the MRCA of all individuals of the initial population?
        :param current_generation: string
        :return: boolean
        """
        
        counter = 0;
        for i in current_generation:
            if(i != "-"):
                counter += 1

        if(counter == 1):
            return True
        else:
            return False

    def make_initial_generation(self):
        """
        makes the initial generation, one letter per individual
        :return: string
        """

        return string.ascii_lowercase[0:self.size]

    def simulate(self):

        # generate initial population
        generation = 0
        waiting_time_overlap = 0
        k = self.size
        current_generation = self.make_initial_generation()
        
        print(str(generation), end="")
        print("\t", end="")
        print(current_generation)

        # run the simulator, print out each generation 
        # until the MRCA of all existing individuals is reached
        while(not self.is_MRCA_of_all_found(current_generation)):
            generation -= 1
            [current_generation, k, waiting_time_overlap] = self.simulate_previous_generation(current_generation, k, waiting_time_overlap)

            print(str(generation), end="")
            print("\t", end="")
            print(current_generation)

        return self.time

def main():
    try:

        args = create_parser()

        # simulate "repeats" times for each 2N = "sizes"
        simulation_sizes = args.sizes
        simulation_repeats = args.repeats
        console_output = args.consoleOutput

        if(not console_output):
            path = './Robin_Bonkass_Sara_Kemmler_generations.txt'
            sys.stdout = open(path, 'w')

        # collect time to find MRCA for every iteration
        # initialize the array used for this
        times = {}
        for size in simulation_sizes:
            times[str(size)] = []


        for size in simulation_sizes:

            print("===========================")

            for i in range(simulation_repeats):

                print(i+1, end="")
                print(". simulation for N2 = ", end="")
                print(size, end="")
                print(":")
                print("---------------------------")

                pop_gen_simulator = PopGenSimulator(size)
                times[str(size)].append(pop_gen_simulator.simulate())

                print("===========================")

            print("\n")

        # get the values for the y axis for the plot
        values = []
        for size in times:
            sum = 0
            for i in times[size]:
                sum += i
            avarage = sum/simulation_repeats
            values.append(avarage)

        # get the theory values
        values_theory = []
        for size in times:
            values_theory.append(2 * (1 - (1 / int(size))) * int(size))
        
        # get the names for the x-Axis for the plot
        names = []
        for size in times:
            name = str(size)
            names.append(name)
        x_axis = np.arange(len(names))

        plt.bar(x_axis - 0.1, values,        0.2, label="avarage values")
        plt.bar(x_axis + 0.1, values_theory, 0.2, label="theory values")

        plt.xticks(x_axis, names)
        plt.xlabel("Amount of Genes 2N")
        plt.ylabel("time to find MRCA of all genes")
        plt.title("time for all initial individuals to find their MRCA\n depending on the population size")
        plt.legend()
        plt.savefig('Robin_Bonkass_Sara_Kemmler_time_to_MRCA.pdf')

    except NameError as err:
        print(err)
    except:
        print('python Robin_Bonkass_Sara_Kemmler_A6.py -fileOutput -repeats 3  -sizes 4 6 10')


if __name__ == "__main__":
    main()
