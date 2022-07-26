'''
Assignment 12
Robin Bonkaß
Sara Kemmler
'''

import argparse

def create_parser():
    '''
    Creates a argument parser to facilitate file reading.
    '''

    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-ftrue', '--fileTrue',
                   help="Fasta file containing true secondary structure")
    p.add_argument('-fphd', '--filePhd',
                   help="Fasta file containing prediction of PHD")
    p.add_argument('-fproteus2', '--fileProteus2',
                   help="Fasta file containing prediction of PHD")
    
    return(p.parse_args())


def parseFileTrue(file):
    '''
    input: file (relies on the format given in the trueSecStructure.txt file)
    returns:  string containing the secondary structure
    '''
    secStructure = ""
    with open(file) as file_content:
        for line in file_content:
            if line.startswith("STR"):
                secStructure += line[4:].strip()

    return secStructure

def parseSecStructure(file):
    '''
    input: txt file containing the secondary structure (with line breaks)
    returns:  string containing the secondary structure
    '''
    secStructure = ""
    with open(file) as file_content:
        for line in file_content:
            secStructure += line.strip()

    return secStructure


def toDSSP(secStructure):
    '''
    reduces dssp classification to the symbols {L, E, H}
    '''
    reducedSecStructure = ""
    for i in secStructure:
        if(i in ["C", "S", "T"]):
            reducedSecStructure += "L"
        elif(i in ["E", "B"]):
            reducedSecStructure += "E"
        elif(i in ["H", "G", "I"]):
            reducedSecStructure += "H"
        elif(i in ["L"]):
            reducedSecStructure += "L"
        else:
            raise Exception("Secondary structure containing unknown symbol " + i)

    return reducedSecStructure


def getQ3(secStructure1, secStructure2):
    '''
    returns the Q3 value for the two given secondary structures
    '''
    if(len(secStructure1) == len(secStructure2)):
        equalCount = 0
        for i in range(len(secStructure1)):
            if(secStructure1[i] == secStructure2[i]):
                equalCount += 1
        Q3 = (equalCount/len(secStructure1))*100
        return Q3
    else:
        raise Exception("The two secondary structures don't have the same length")



def main():
    '''
    The main function
    '''

    fileTrue        = args.fileTrue
    filePhd         = args.filePhd
    fileProteus         = args.fileProteus2

    secStructureTrue     = toDSSP(parseFileTrue(fileTrue))
    secStructurePhd      = toDSSP(parseSecStructure(filePhd))
    secStructureProteus  = toDSSP(parseSecStructure(fileProteus))

    print("\nQ3 of the true secondary structure and the PHD prediction:")
    print(str(getQ3(secStructureTrue, secStructurePhd)) + "%")
    print("\nQ3 of the true secondary structure and the Proteus2 prediction:")
    print(str(getQ3(secStructureTrue, secStructureProteus)) + "%")


if __name__ == "__main__":
    # try:
        args = create_parser()
        
        # accessing the arguments of argparse
        print("True Secondary structure: " + args.fileTrue)
        print("PHD prediction: " + args.filePhd)
        print("Proteus2 prediction: " + args.fileProteus2)

        main()
    # except:
    #     print('Try:  python3 Sara_Kemmler_Robin_Bonkass_A12.py -ftrue trueSecStructure.txt -fphd secStructure_PHD.txt -fproteus2 secStructure_proteus.txt')