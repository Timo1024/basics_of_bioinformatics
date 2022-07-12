from numpy import number
from gff import gffParser
from itertools import chain
import sys

def main():

    out1 = gffParser(input_file)
    out2 = gffParser(input_file2)


    # print(out.getCDS("NC_002516.2", 1))
    # print(out.getGenes("NC_002516.2"))

    # get all starts and ends
    # start_end1 = []
    numbers1 = []
    for i in out1.getCDS("NC_002516.2", "gene-PA0001"):
        # start_end1.append([i['start'], i['end']])
        # numbers1 += list(range(i['start'], i['end']+1))
        numbers1 = chain(numbers1, range(i['start'], i['end']+1))

    # start_end2 = []
    numbers2 = []
    for i in out2.getCDS("NC_002516.2", "gene-PA0001"):
        # start_end2.append([i['start'], i['end']])
        # numbers2 += list(range(i['start'], i['end']+1))
        numbers2 = chain(numbers2, range(i['start'], i['end']+1))

    # print(numbers1)

    # z = [i for i in numbers1 if i in numbers2]
    # print(z)

    set_numbers1 = set(numbers1)
    print(set_numbers1.intersection(numbers2))


    # numbers1 = range(5,10)
    # numbers1 = chain(numbers1, range(10,15))
    # numbers2 = chain(range(1,7), range(12,15))
    # xs = set(numbers1)
    # print(xs.intersection(numbers2))

if __name__ == "__main__":
    # try:
        input_file = sys.argv[1]
        input_file2 = sys.argv[2]
        print(input_file)
        print(input_file2)
        main()
    # except:
    #     print('Try:  python3 Robin_Bonkass_Sara_Kemmler_A10.py PAO1_annotation.gff')