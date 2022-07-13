
def getLines(file):
    '''
    gets all important lines from the .gff file
    '''
    lines = []
    with open(file) as file_content:
        for line in file_content:
            if(not line.startswith("#") and len(line) > 1):
                lines.append(line)
    return lines

def getCDSStartEndFromLine(line):
    '''
    gets the start and end of an cds from one line if it exists
    otherwise returns False
    '''
    line_array = line.strip().split("\t")
    if(len(line_array) >= 3):
        index = 0
        for i in line_array:
            if i == "CDS":
                return [int(line_array[index+1]), int(line_array[index+2])]
            index += 1
        return False
    else:
        return False

def getAllStartEnd(file):
    '''
    gets all start end end points of all existing cds
    and returns them as an array
    '''
    all_start_end = []
    lines = getLines(file)
    for i in lines:
        if getCDSStartEndFromLine(i) != False:
            all_start_end.append(getCDSStartEndFromLine(i))
    return all_start_end