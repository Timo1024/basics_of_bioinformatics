from cmath import exp
import math
import string
import sys
import random  # needed in simulate_previous_generation

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

    # print(first_letter_index)
    # print(second_letter_index)

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
        if(size < 0 or size > 26):
            print("Parameter size has to be an int between 0 and 26")
            self.size = 0  # population size
        else:
            self.size = size  # population size

    # def get_random_pair(current_generation_array):
    #     """
    #     returns array with a "-" at one position instead of a letter
    #     -> this means that two letters are merged two one
    #     """

    #     length = len(current_generation_array)

    #     # counts all left letters
    #     letters_count = 0
    #     for i in range(length):
    #         if(current_generation_array[i] != "-"):
    #             letters_count += 1

    #     # get index for two letters which will be merged
    #     first_letter_index = math.floor(random.random() * letters_count)
    #     second_letter_index = math.floor(random.random() * letters_count)
    #     if(first_letter_index == second_letter_index):
    #         second_letter_index += 1
    #         second_letter_index = second_letter_index % letters_count

    #     # get first letter
    #     first = 0
    #     for i in range(length):
    #         if(first == first_letter_index):
    #             first_letter = current_generation_array[i]
    #         else:
    #             if(current_generation_array != "-"):
    #                 first += 1

    #     # get second letter
    #     second = 0
    #     for i in range(length):
    #         if(second == second_letter_index):
    #             second_letter = current_generation_array[i]
    #         else:
    #             if(current_generation_array != "-"):
    #                 second += 1

        # # get lexicographic lower letter -> save this as first_letter and other as second_letter
        # if(first_letter > second_letter):
        #     first_letter_temp = first_letter
        #     first_letter = second_letter
        #     second_letter = first_letter_temp

        # # get position of second letter and make it a "-"
        # merged_letter_index = current_generation_array.index(second_letter)
        # current_generation_array[merged_letter_index] = "-"

        # return current_generation_array

    # TODO needs k???
    def simulate_previous_generation(self, current_generation, k):
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
        # TODO PLEASE IMPLEMENT

        # make array from string
        current_generation_array=[]
        current_generation_array[:0]=current_generation
        # current_generation_array = list(current_generation)

        waiting_time = 0

        while(waiting_time<=1 and k >= 2):
            random_uniform = random.uniform(0,1)
            waiting_time += ((0-(math.comb(k,2) ** (-1))) * math.log(random_uniform)) * k
            if(waiting_time <= 1):
                # print(waiting_time)
                # print(k)
                current_generation_array = get_random_pair(current_generation_array)
                k -= 1

        previous_generation = "".join(current_generation_array)

        return [previous_generation, k]


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
        # FIXME deleted size as parameter
        """
        makes the initial generation, one letter per individual
        :param size: int
        :return: string
        """

        return string.ascii_lowercase[0:self.size]


def main():
    # try:
        size = int(sys.argv[1])

        # create instance of class PopGenSimulator
        pop_gen_simulator = PopGenSimulator(size)

        # generate initial population
        generation = 0
        k = size
        current_generation = pop_gen_simulator.make_initial_generation()
        print(str(generation), end="")
        print("\t", end="")
        print(current_generation)

        # run the simulator, print out each generation 
        # until the MRCA of all existing individuals is reached
        while(not pop_gen_simulator.is_MRCA_of_all_found(current_generation)):
        # for i in range(size-1):
            generation -= 1
            [current_generation, k] = pop_gen_simulator.simulate_previous_generation(current_generation, k)
            # k -= 1
            print(str(generation), end="")
            print("\t", end="")
            print(current_generation)

    # except:
    #     print('run name1_name2_PopGenSimulator.py 4')


if __name__ == "__main__":
    main()
