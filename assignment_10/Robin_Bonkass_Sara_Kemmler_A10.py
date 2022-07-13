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
    cds1 = out1.getCDS("NC_002516.2", "gene-PA0001")
    cds2 = out2.getCDS("NC_002516.2", "gene-PA0001")

    # gets the max length of sequence
    last_cds1_end = cds1[len(cds1)-1]['end']
    last_cds2_end = cds2[len(cds2)-1]['end']
    last_cds_end_all = max(last_cds1_end, last_cds2_end)
    
    # gets the ranges for cds in ground truth
    ground_truth = []
    for i in cds1:
        # start_end1.append([i['start'], i['end']])
        ground_truth += list(range(i['start'], i['end']+1))
        # ground_truth = chain(ground_truth, range(i['start'], i['end']+1))

    # gets the ranges for cds in prediction
    prediction = []
    for i in cds2:
        # start_end2.append([i['start'], i['end']])
        prediction += list(range(i['start'], i['end']+1))
        # prediction = chain(prediction, range(i['start'], i['end']+1))

    
    ground_truth_bool = [False] * (last_cds_end_all+1)
    prediction_bool = [False] * (last_cds_end_all+1)

    for i in ground_truth:
        ground_truth_bool[i] = True
    for i in prediction:
        prediction_bool[i] = True

    tp = 0
    tn = 0
    fp = 0
    fn = 0

    for i in range(last_cds_end_all+1):
        if(ground_truth_bool[i] and prediction_bool[i]):
            tp += 1
        elif ground_truth_bool[i]:
            fn += 1
        elif prediction_bool[i]:
            fp += 1
        else:
            tn += 1

    print("\nAfter the first " + str(i) + " nucleotides, the values are:")
    print("tp = " + str(tp))
    print("tn = " + str(tn))
    print("fp = " + str(fp))
    print("fn = " + str(fn))
    if(tp == 0 and fn == 0):
        print("sensitivity = \t" + "undefined")
    else:
        print("sensitivity = \t" + str(tp/(tp+fn)))
    if(fp == 0 and tn == 0):
        print("specificity = \t" + "undefined")
    else:
        print("specificity = \t" + str(tn/(fp+tn)))
    print("accuracy = \t" + str((tp+tn)/(i+1)))


if __name__ == "__main__":
    # try:
        input_file = sys.argv[1]
        input_file2 = sys.argv[2]
        print(input_file)
        print(input_file2)
        main()
    # except:
    #     print('Try:  python3 Robin_Bonkass_Sara_Kemmler_A10.py PAO1_annotation.gff Galaxy2-[Prokka_on_data_1__gff].gff3')