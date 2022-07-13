# TODO comments

def getLines(file):
    lines = []
    with open(file) as file_content:
        for line in file_content:
            if(not line.startswith("#") and len(line) > 1):
                lines.append(line)
    return lines

def getCDSStartEndFromLine(line):
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
    all_start_end = []
    lines = getLines(file)
    for i in lines:
        # print(i)
        # print(getCDSStartEndFromLine(i))
        if getCDSStartEndFromLine(i) != False:
            all_start_end.append(getCDSStartEndFromLine(i))
    return all_start_end