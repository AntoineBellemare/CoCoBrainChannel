import numpy as np
import math
from fractions import Fraction
import itertools
from biotuner_utils import *

################################# PEAK RATIOS AND HARMONICS ####################################################
#This function calculates the 15 ratios (with the possibility to bound them between 1 and 2) derived from the 6 peaks
def compute_peak_ratios(peaks, rebound = 1, octave = 2):
    ratios = []
    peak_ratios_rebound = []
    for p1 in peaks:
        for p2 in peaks:
            ratio_temp = p2/p1
            if ratio_temp == 1:
                ratio_temp = None
            elif ratio_temp < 1:
                ratio_temp = None

            ratios.append(ratio_temp)

        peak_ratios = np.array(ratios)
        peak_ratios = [i for i in peak_ratios if i]
        peak_ratios = list(set(peak_ratios))
    if rebound == 1:
        for peak in peak_ratios:
            while peak > octave:
                peak = peak/octave
            peak_ratios_rebound.append(peak)
    peak_ratios_rebound = np.array(peak_ratios_rebound)
    peak_ratios_rebound = [i for i in peak_ratios_rebound if i]
    peak_ratios = sorted(peak_ratios)
    peak_ratios_rebound = sorted(list(set(peak_ratios_rebound)))
    return peak_ratios, peak_ratios_rebound

    #This function takes a list of frequency peaks as input and computes the desired number of harmonics
#with the formula: x + 2x + 3x ... + nx
def EEG_harmonics_mult(peaks, n_harmonics, n_oct_up = 0):
    n_harmonics = n_harmonics + 2
    multi_harmonics = []
    multi_harmonics_rebound = []
    for p in peaks:
        multi_harmonics_r = []
        multi_harm_temp = []
        harmonics = []
        p = p * (2**n_oct_up)
        i = 1
        harm_temp = p
        while i < n_harmonics:
            harm_temp = p * i
            harmonics.append(harm_temp)
            i+=1
        multi_harmonics.append(harmonics)
    multi_harmonics = np.array(multi_harmonics)

    return multi_harmonics

#This function takes a list of frequency peaks as input and computes the desired number of harmonics
#with the formula: x + x/2 + x/3 ... + x/n
def EEG_harmonics_div(peaks, n_harmonics, n_oct_up = 0):
    n_harmonics = n_harmonics + 2
    multi_harmonics = []
    multi_harmonics_sub = []
    for p in peaks:

        harmonics = []
        harmonics_sub = []
        p = p * (2**n_oct_up)
        i = 2
        harm_temp = p
        harm_temp_sub = p
        while i < n_harmonics:
            harm_temp = harm_temp + (p/i)
            harm_temp_sub = abs(harm_temp_sub - (p/i))
            harmonics.append(harm_temp)
            harmonics_sub.append(harm_temp_sub)
            i+=1
        multi_harmonics.append(harmonics)
        multi_harmonics_sub.append(harmonics_sub)
    multi_harmonics = np.array(multi_harmonics)
    multi_harmonics_bounded = multi_harmonics.copy()
    multi_harmonics_sub = np.array(multi_harmonics_sub)
    multi_harmonics_sub_bounded = multi_harmonics_sub.copy()
    #Rebound the result between 1 and 2
    for i in range(len(multi_harmonics_bounded)):
        for j in range(len(multi_harmonics_bounded[0])):
            multi_harmonics_bounded[i][j] = rebound(multi_harmonics_bounded[i][j])
            multi_harmonics_sub_bounded[i][j] = rebound(multi_harmonics_sub_bounded[i][j])
    return multi_harmonics, multi_harmonics_bounded, multi_harmonics_sub, multi_harmonics_sub_bounded

#This function computes harmonics of a list of peaks and compares the lists of harmonics pairwise to find fitting
#between the harmonic series
def harmonic_fit(peaks, n_harm = 10, bounds = 1, function = 'mult'):
    from itertools import combinations
    peak_bands = []
    for i in range(len(peaks)):
        peak_bands.append(i)
    if function == 'mult':
        multi_harmonics = EEG_harmonics_mult(peaks, n_harm)
    elif function == 'div':
        multi_harmonics, x, y, z = EEG_harmonics_div(peaks, n_harm)
    #print(multi_harmonics)
    list_peaks = list(combinations(peak_bands,2))
    #print(list_peaks)
    harm_temp = []
    for i in range(len(list_peaks)):
        harms, b, c = compareLists(multi_harmonics[list_peaks[i][0]], multi_harmonics[list_peaks[i][1]], bounds)
        harm_temp.append(harms)
    harm_fit = np.array(harm_temp).squeeze()

    if len(peak_bands) > 2:
        harm_fit = list(itertools.chain.from_iterable(harm_fit))
        harm_fit = [round(num, 3) for num in harm_fit]
        harm_fit = list(dict.fromkeys(harm_fit))
    return harm_fit



    ####################################   N-TET   ##############################################################


    #Oct_subdiv
#Argument 1 : a ratio in the form of a float or a fraction
#Argument 2 : bounds between which the octave should fall
#Argument 3 : value of the octave
#Argument 4 : number of octave subdivisions

def oct_subdiv(ratio,bounds,octave,n):
    Octdiv, Octvalue, i = [], [], 1
    ratios = []
    while len(Octdiv) < n:
        ratio_mult = (ratio**i)
        while ratio_mult > octave:
            ratio_mult = ratio_mult/octave

        rescale_ratio = ratio_mult - round(ratio_mult)
        ratios.append(ratio_mult)
        i+=1
        if -bounds < rescale_ratio < bounds:
            Octdiv.append(i-1)
            Octvalue.append(ratio_mult)
        else:
            continue
    return Octdiv, Octvalue, ratios



def compare_oct_div(Octdiv = 53, Octdiv2 = 12, bounds = 0.01, octave = 2):
    ListOctdiv = []
    ListOctdiv2 = []
    OctdivSum = 1
    OctdivSum2 = 1
    i = 1
    i2 = 1
    Matching_harmonics = []
    #HARMONIC_RATIOS = [1, 1.0595, 1.1225, 1.1892, 1.2599, 1.3348, 1.4142, 1.4983, 1.5874, 1.6818, 1.7818, 1.8897]
    while OctdivSum < octave:
        OctdivSum =(nth_root(octave, Octdiv))**i
        i+=1
        ListOctdiv.append(OctdivSum)
    #print(ListOctdiv)

    while OctdivSum2 < octave:
        OctdivSum2 =(nth_root(octave, Octdiv2))**i2
        i2+=1
        ListOctdiv2.append(OctdivSum2)
    #print(ListOctdiv2)
    for i, n in enumerate(ListOctdiv):
        for j, harm in enumerate(ListOctdiv2):
            if harm-bounds < n < harm+bounds:
                Matching_harmonics.append([n, i+1, harm, j+1])
    Matching_harmonics = np.array(Matching_harmonics)
    return Matching_harmonics


##Here we use the most consonant peaks ratios as input of oct_subdiv function. Each consonant ratio
##leads to a list of possible octave subdivisions. These lists are compared and optimal octave subdivisions are
##determined.

#Output1: octave subdivisions
#Output2: ratios that led to Output1

def multi_oct_subdiv (peaks, max_sub, octave_limit):
    import itertools
    from collections import Counter
    a, b = compute_consonance(peaks, 0.01)
    c, ratios = comp_consonant_integers_ratios(b)

    list_oct_div = []
    for i in range(len(ratios)):
        list_temp, no, no2 = oct_subdiv(ratios[i], octave_limit, 2, 30)
        list_oct_div.append(list_temp)


    counts = Counter(list(itertools.chain(*list_oct_div)))
    oct_div_temp = []
    for k, v in counts.items():
        if v > 1:
            oct_div_temp.append(k)
    oct_div_temp = np.sort(oct_div_temp)
    oct_div_final = []
    for i in range(len(oct_div_temp)):
        if oct_div_temp[i] < max_sub:
            oct_div_final.append(oct_div_temp[i])
    return oct_div_final, ratios




    #######################################   CONSONANCE   ###########################################################



def compute_consonance (peaks, limit):
    from fractions import Fraction
    consonance_ = []
    peaks2keep = []
    peaks_consonance = []
    for p1 in peaks:
        for p2 in peaks:
            peaks2keep_temp = []
            p2x = p2
            p1x = p1
            if p1x > p2x:
                while p1x > p2x:
                    p1x = p1x/2
            if p1x < p2x:
                while p2x > p1x:
                    p2x = p2x/2

            cons_temp = Fraction(p2x/p1x).limit_denominator(1000)
            #print(cons_temp)
            cons_temp = (cons_temp.numerator + cons_temp.denominator)/(cons_temp.numerator * cons_temp.denominator)
            #print(cons_temp)
            #if limit > 0:
            if cons_temp >= 1 :
                cons_temp = None
                p2x = None
                p1x = None

            elif cons_temp < limit:
                cons_temp = None
                p2x = None
                p1x = None
            if p2x != None:
                peaks2keep_temp.extend([p2, p1])
            consonance_.append(cons_temp)
            peaks2keep.append(peaks2keep_temp)
        peaks_consonance = np.array(peaks2keep)
        peaks_consonance = [x for x in peaks_consonance if x]
        consonance = np.array(consonance_)
        consonance = [i for i in consonance if i]

        #consonance = list(set(consonance))
    return consonance, peaks_consonance

## Function that keeps the frequencies that are the most consonant with others
##Takes pairs of frequencies that are consonant (output of the 'compute consonance' function)
def multi_consonance(pairs, n_freqs = 5):
    freqs_dup = list(itertools.chain(*pairs))
    pairs_temp = list(itertools.chain.from_iterable(pairs))
    freqs_nodup = list(dict.fromkeys(pairs_temp))
    f_count = []
    for f in freqs_nodup:
        f_count.append(freqs_dup.count(f))
    freqs_related = [x for _,x in sorted(zip(f_count,freqs_nodup))][-(n_freqs):][::-1]
    return freqs_related


# Function that computes integer ratios from peaks with higher consonance
# Needs at least two pairs of values
def comp_consonant_integers_ratios (peaks_consonance):
    from fractions import Fraction
    cons_integer = []

    for i in range(len(peaks_consonance)):
        #for j in range(len(peaks_consonance[0])):
        a = peaks_consonance[i][0]
        b = peaks_consonance[i][1]
        if a > b:
            while a > (b):
                a = a/2
        if a < b:
            while b > (a):
                b = b/2
        #c = a/b
        cons_temp = Fraction(a/b).limit_denominator(1000)
        num = cons_temp.numerator
        denom = cons_temp.denominator
        frac_list = [num, denom]
        frac = tuple(frac_list)
        #print(cons_temp)
        cons_integer.append(frac)
    consonant_integers = np.array(cons_integer).squeeze()
    cons_ratios = []
    for j in range(len(consonant_integers)):
        cons_ratios.append((consonant_integers[j][0])/(consonant_integers[j][1]))
    consonant_ratios = np.array(cons_ratios)
    cons_rat = []
    for i in consonant_ratios:
        if i not in cons_rat:
            cons_rat.append(i)
    try:
        cons_rat.remove(1.0)
    except (ValueError, TypeError, NameError):
        pass
        #cons_int = []
        #for i in consonant_integers:
        #    if i not in cons_int:
        #        cons_int.append(i)
    return consonant_integers, cons_rat


##Require a function that computes
