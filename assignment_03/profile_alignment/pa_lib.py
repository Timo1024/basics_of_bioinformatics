import argparse
from pa_classes import *
import numpy as np

'''
Here are all functions which are used in the main file
'''

def create_parser():
    '''
    Creates a argument parser to facilitate file reading.

    '''
    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-f1', '--file-one',
                   help="File for the the programming task")

    return(p.parse_args())


def read_arguments():
    '''
    reads the match/mismatch and gap penalty from the command line
    '''

    input_1 = False
    input_2 = False
    input_3 = False

    match = 0
    mismatch = 0
    gap = 0

    # gets input for matching score
    while(not input_1):
        match = input("Type the score for a matching pair:")
        try:
            match = int(match)
            input_1 = True
        except:
            print("Wrong input")

    # gets input for mismatching score
    while(not input_2):
        mismatch = input("Type the score for a mismatch:")
        try:
            mismatch = int(mismatch)
            input_2 = True
        except:
            print("Wrong input")

    # gets input for gap penalty
    while(not input_3):
        gap = input("Type the score for gap penalty:")
        try:
            gap = int(gap)
            input_3 = True
        except:
            print("Wrong input")

    return [match, mismatch, gap]


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

    # concats all line which represent one sequence
    sequence_array_concat = []
    for sequence in sequence_array:
        seq = ""
        for line in sequence:
            seq += line
        sequence_array_concat.append(seq)

    return [heading_array, sequence_array_concat]


def initialize_paList(sequences):
    '''
    makes profile alignment (pa) objects from each sequence which is given
    so every pa contains of one unchanged sequence
    gives the initial profile alignment list array back
    (paList)
    '''

    paList = []
    for seq in sequences:
        paList.append(profile_alignment([seq]))
    return paList


def calc_scross(paList, s, d):
    '''
    calculates the maximum S_cross score of a given paList
    gives back the paList without the two pa's which made the
    highest cross score
    and gives back the two pa's and the sequences of them which
    made this score
    returns [paList, [pa1, sequence_of_pa1], [pa2, sequence_of_pa2]]
    '''

    # makes the crossproduct paList x paList
    # and compares every tuple for the S_cross score
    max_pas = []
    max_seqs = []
    max_scores = []
    for pa1 in paList:
        for pa2 in paList:
            if(not pa1.getID() == pa2.getID()):

                # calculate S_cross for the tuple (pa1, pa2)
                seqs = []
                scores = []
                for seq1 in pa1.seqs:
                    for seq2 in pa2.seqs:
                        new_alignment = pair_alignment(seq1, seq2, s, d)
                        score = new_alignment.getScore()
                        seqs.append([seq1, seq2])
                        scores.append(score)
                max_score_index = np.argmax(scores)
                max_pas.append([pa1, pa2])
                max_seqs.append(seqs[max_score_index])
                max_scores.append(scores[max_score_index])
    
    # look which score was the highest
    max_score_index = np.argmax(max_scores)
    newpaList = []

    first_occurance_1 = False
    first_occurance_2 = False

    for pa in paList:
        if((pa.getID() != max_pas[max_score_index][0].getID() or first_occurance_1) and (pa.getID() != max_pas[max_score_index][1].getID() or first_occurance_2)):
            newpaList.append(pa)

        # checks if one of the used pas was detected and therefore not added to the newpaList
        # sets the first occurance to True, so that if there is another equal sequence
        # it adds it to the newpaList b/c we just want to delete the one we used for the max
        # S_cross value
        if(pa.getID() == max_pas[max_score_index][0].getID()):
            first_occurance_1 = True
        if(pa.getID() == max_pas[max_score_index][1].getID()):
            first_occurance_2 = True

    return [newpaList, [max_pas[max_score_index][0], max_seqs[max_score_index][0]], [max_pas[max_score_index][1], max_seqs[max_score_index][1]]]

            
