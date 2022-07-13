from gff_parser import getAllStartEnd
import sys

def getValues(ground_truth_file, prediction_file):

    # get all starts and ends
    cds1 = getAllStartEnd(ground_truth_file)
    cds2 = getAllStartEnd(prediction_file)

    # gets the max length of sequence
    last_cds1_end = cds1[len(cds1)-1][1]
    last_cds2_end = cds2[len(cds2)-1][1]
    last_cds_end_all = max(last_cds1_end, last_cds2_end)
    
    # gets the ranges for cds in ground truth
    ground_truth = []
    for i in cds1:
        ground_truth += list(range(i[0], i[1]+1))

    # gets the ranges for cds in prediction
    prediction = []
    for i in cds2:
        prediction += list(range(i[0], i[1]+1))

    # boolean arrays. True if index is in cds
    ground_truth_bool = [False] * (last_cds_end_all+1)
    prediction_bool = [False] * (last_cds_end_all+1)

    for i in ground_truth:
        ground_truth_bool[i] = True
    for i in prediction:
        prediction_bool[i] = True

    # calculating tp, tn, fp, fn
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

    sensitivity = tp/(tp+fn) or 0
    specificity = tn/(fp+tn) or 0
    accuracy    = (tp+tn)/(i+1)

    return [tp, fn, fp, tn, sensitivity, specificity, accuracy]

def printResults(results):
    [tp, fn, fp, tn, sensitivity, specificity, accuracy] = results

    print("\nThe calculated values are:")
    print("tp = " + str(tp))
    print("tn = " + str(tn))
    print("fp = " + str(fp))
    print("fn = " + str(fn))
    print("sensitivity = \t" + str(sensitivity))
    print("specificity = \t" + str(specificity))
    print("accuracy = \t" + str(accuracy))

def main():

    # for file1 and file2
    values = getValues(input_file, input_file2)

    print("\nValues for file " + input_file + " as ground truth and " + input_file2 + " as prediction")
    printResults(values)

    # for file1 and file3
    values = getValues(input_file, input_file3)

    print("\n\nValues for file " + input_file + " as ground truth and " + input_file3 + " as prediction")
    printResults(values)

    # # get all starts and ends
    # cds1 = getAllStartEnd(input_file)
    # cds2 = getAllStartEnd(input_file2)

    # # gets the max length of sequence
    # last_cds1_end = cds1[len(cds1)-1][1]
    # last_cds2_end = cds2[len(cds2)-1][1]
    # last_cds_end_all = max(last_cds1_end, last_cds2_end)
    
    # # gets the ranges for cds in ground truth
    # ground_truth = []
    # for i in cds1:
    #     ground_truth += list(range(i[0], i[1]+1))

    # # gets the ranges for cds in prediction
    # prediction = []
    # for i in cds2:
    #     prediction += list(range(i[0], i[1]+1))

    # # boolean arrays. True if index is in cds
    # ground_truth_bool = [False] * (last_cds_end_all+1)
    # prediction_bool = [False] * (last_cds_end_all+1)

    # for i in ground_truth:
    #     ground_truth_bool[i] = True
    # for i in prediction:
    #     prediction_bool[i] = True

    # # calculating tp, tn, fp, fn
    # tp = 0
    # tn = 0
    # fp = 0
    # fn = 0

    # for i in range(last_cds_end_all+1):
    #     if(ground_truth_bool[i] and prediction_bool[i]):
    #         tp += 1
    #     elif ground_truth_bool[i]:
    #         fn += 1
    #     elif prediction_bool[i]:
    #         fp += 1
    #     else:
    #         tn += 1

    # print("\nAfter the first " + str(i) + " nucleotides, the values are:")
    # print("tp = " + str(tp))
    # print("tn = " + str(tn))
    # print("fp = " + str(fp))
    # print("fn = " + str(fn))
    # if(tp == 0 and fn == 0):
    #     print("sensitivity = \t" + "undefined")
    # else:
    #     print("sensitivity = \t" + str(tp/(tp+fn)))
    # if(fp == 0 and tn == 0):
    #     print("specificity = \t" + "undefined")
    # else:
    #     print("specificity = \t" + str(tn/(fp+tn)))
    # print("accuracy = \t" + str((tp+tn)/(i+1)))


if __name__ == "__main__":
    # try:
        input_file = sys.argv[1]
        input_file2 = sys.argv[2]
        input_file3 = sys.argv[3]
        main()
    # except:
    #     print('Try:  python3 Robin_Bonkass_Sara_Kemmler_A10.py PAO1_annotation.gff Galaxy2-[Prokka_on_data_1__gff].gff3')