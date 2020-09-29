import numpy as np


def rebound(x, low = 1, high = 2, octave = 2):
    while x > high:
        x = x/octave
    while x < low:
        x = x*octave
    return x

def nth_root (num, root):
    answer = num**(1/root)
    return answer


#Function that compares lists (i.e. peak harmonics)
def compareLists(list1, list2, bounds):
    matching = []
    matching_pos = []
    for i, l1 in enumerate(list1):
        for j, l2 in enumerate(list2):
            if l2-bounds < l1 < l2+bounds:
                matching.append((l1+l2)/2)
                matching_pos.append([(l1+l2)/2, i+1, j+1])
    matching = np.array(matching)
    matching_pos = np.array(matching_pos)
    ratios_temp = []
    for i in range(len(matching_pos)):
        if matching_pos[i][1]>matching_pos[i][2]:
            ratios_temp.append(matching_pos[i][1]/matching_pos[i][2])
            #print(matching_pos[i][0])
        else:
            ratios_temp.append(matching_pos[i][2]/matching_pos[i][1])
    matching_pos_ratios = np.array(ratios_temp)
    return matching, matching_pos, matching_pos_ratios

def create_SCL(scale, fname):
    #Output SCL files
    from pytuning.scales.pythagorean import create_pythagorean_scale
    from pytuning.tuning_tables import create_scala_tuning
    #scale = create_pythagorean_scale()
    table = create_scala_tuning(scale,fname)
    outF = open(fname+'.scl', "w")
    outF.writelines(table)
    outF.close()
    return
