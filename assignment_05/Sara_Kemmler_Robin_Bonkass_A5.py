'''
Assignment 5
Robin Bonkaß
Sara Kemmler
'''

import argparse
import copy
from math import sqrt

import numpy as np

def print_matrix(matrix, full = False):
    '''
    prints a matrix in command line
    used to debug
    '''
    
    for i in range(len(matrix)):
        if(not full):
            for j in range(len(matrix[i])):
                print(str(round(matrix[i][j], 2)), "\t", end = '')
            print("")
        else:
            for j in range(len(matrix[i])):
                if(j == 0 and i == 0):
                    print("\t", end = '')
                else:
                    print(str(round(matrix[i][j], 2)), "\t", end = '')
            print("")

def create_parser():
    '''
    Creates a argument parser to facilitate file reading.
    '''

    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-f1', '--file-one',
                   help="Original Distances")
    p.add_argument('-f2', '--file-two',
                   help="Distances tree 1")
    p.add_argument('-f3', '--file-three',
                   help="Distances tree 2")

    return(p.parse_args())

def create_matrix(file):
    matrix = []
    file_length = 0

    with open(file) as file_content:
        for line in file_content:
            file_length += 1

    with open(file) as file_content:
        line_count = 0
        for line in file_content:
            line_array = line.split("\n")
            # FIXME doesn't work when seperated by spaces instead of tabs
            line_array = line_array[0].split("\t")[1:line_count+1]
            line_array = [float(x) for x in line_array]
            line_array.append(0)
            matrix.append(line_array)

            for i in range(file_length - line_count - 1):
                matrix[line_count].append(-1)

            line_count += 1

        # complete other half of matrix
        length = len(matrix)
        for i in range(length):
            for j in range(len(matrix[i])):
                matrix[i][j] = matrix[j][i]

    return matrix

# TODO ist das was mit D Strich und T Strich gemeint ist?
# FIXME sollen hier nur die Zahlen mit einberechnet werden, die in dem file stehen?
# Berechnet zur Zeit durchschnitt über alle Zahlen, auch die 0 in der Diagonalen
def avarage_distances(matrix):
    '''
    returns the avarage value of all entries of the matrix
    '''

    n = len(matrix)
    sum = 0

    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if(j > i):
                sum += matrix[i][j]

    avarage = sum * (2/(n*(n-1)))

    return avarage


def ccc(D, T):
    '''
    calculates the cophenetic correlation coefficient of
    the distance matrices D and T
    '''

    avarage_D = avarage_distances(D)
    avarage_T = avarage_distances(T)

    # calculating the numerator
    numerator = 0
    for i in range(len(D)):
        for j in range(len(D[i])):
            if(i<j):
                numerator += ((D[i][j] - avarage_D) * (T[i][j] - avarage_T))

    # calculating the denominator
    denominator = 0

    sum_D = 0
    for i in range(len(D)):
        for j in range(len(D[i])):
            if(i<j):
                sum_D += ((D[i][j] - avarage_D)**2)

    sum_T = 0
    for i in range(len(T)):
        for j in range(len(T[i])):
            if(i<j):
                sum_T += ((T[i][j] - avarage_T)**2)

    denominator = sqrt(sum_D * sum_T)

    return numerator/denominator


def main():
    '''
    The main function should contain all functions that solved the stated tasks. 
    '''

    file_original = args.file_one
    file_tree1    = args.file_two
    file_tree2    = args.file_three

    original = create_matrix(file_original)
    tree1    = create_matrix(file_tree1)
    tree2    = create_matrix(file_tree2)

    print_matrix(original)
    print(avarage_distances(original))
    print("\n")
    print_matrix(tree1)
    print(avarage_distances(tree1))
    print("\n")
    print_matrix(tree2)
    print(avarage_distances(tree2))

    print("")
    print(ccc(original, tree1))

    print("")
    print(ccc(original, tree2))

    
    

if __name__ == "__main__":
    # try:
        args = create_parser()
        
        # accesing the path of the files
        print(args.file_one)
        print(args.file_two)
        print(args.file_three)

        main()
    # except:
    #     print('Try:  python3 Sara_Kemmler_Robin_Bonkass_A5.py -f1 distances_original.dist -f2 distances_tree1.dist -f3 distances_tree2.dist')