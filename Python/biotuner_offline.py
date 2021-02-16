import mne
from mne.datasets.brainstorm import bst_raw
from mne.time_frequency import tfr_morlet, psd_multitaper, psd_welch
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import scipy.io as sio

import pandas as pd
import scipy.io as sio
from scipy.io import savemat, loadmat

#This function is used to compute peak values for each frequency bands on evoked data
def compute_peaks(epochs, condition, chs, FREQ_BANDS, tmin = None, tmax = None, precision = 0.125, sf = 1000):
    import functools
    #epochs are cropped as desire (tmin could be before '0', ex: -1.5, depending on the values used during epoching)
    epochs = epochs.apply_baseline((-1.5, -0.1))
    epochs = epochs.crop(tmin, tmax)
    print(epochs.get_data().shape)

    EOG_chs = ['E1', 'E8', 'E25', 'E32', 'E126', 'E127']
    Unwanted = ['E43', 'E48', 'E49', 'E128', 'E113', 'E120', 'E125', 'E119', 'E129']
    All_chs = epochs.info['ch_names'][0:129]
    EEG_chs = [ele for ele in All_chs if ele not in Unwanted]
    EEG_chs = [ele for ele in EEG_chs if ele not in EOG_chs]
    #Find number of samples
    n_samples = len(epochs.get_data()[0][0])-1
    print(n_samples)
    precision = precision/(sf/n_samples)
    fft_size = int(n_samples/precision)
    print(fft_size)

    evoked = epochs[condition].average(chs)
    FREQs = []
    #This loop iterates for each epoch

    for min, max in FREQ_BANDS:
        #psds, freqs = function(epochs[t], fmin=min, fmax=max, bandwidth = 4, picks = EEG_chs)  #PSDs are calculated with this function, giving power values and corresponding frequency bins as output
        psds, freqs = psd_welch(evoked, fmin=min, fmax=max, n_fft = fft_size)
        psds = 10. * np.log10(psds)   #PSDs values are transformed in log scale to compensate for the 1/f natural slope
        index_max = np.argmax(np.array(psds[13][:]))
        freq = freqs[index_max]
        print(index_max)
        #psds_mean = np.average(psds, axis=1) #Average across bins to obtain one value for the entire frequency range
        FREQs.append(freq)
    FREQs = np.array(FREQs)
    return FREQs

    #This function is used to compute peak values for each frequency bands on evoked data
def compute_peaks_avg(epochs, condition, chs, FREQ_BANDS, dim_reduc = 'avg', tmin = None, tmax = None, precision = 0.125, sf = 1000):
    import functools
    import pandas as pd
    from mne.time_frequency import tfr_morlet, psd_multitaper, psd_welch
    #epochs are cropped as desire (tmin could be before '0', ex: -1.5, depending on the values used during epoching)
    epochs = epochs.apply_baseline((-1.5, -0.1))
    epochs = epochs.crop(tmin, tmax)
    print(epochs.get_data().shape)
    EOG_chs = ['E1', 'E8', 'E25', 'E32', 'E126', 'E127']
    Unwanted = ['E43', 'E48', 'E49', 'E128', 'E113', 'E120', 'E125', 'E119', 'E129']
    All_chs = epochs.info['ch_names'][0:129]
    EEG_chs = [ele for ele in All_chs if ele not in Unwanted]
    EEG_chs = [ele for ele in EEG_chs if ele not in EOG_chs]
    #Find number of samples
    n_samples = len(epochs.get_data()[0][0])-1
    print(n_samples)
    precision = precision/(sf/n_samples)
    fft_size = int(n_samples/precision)
    print(fft_size)
    #Average across trials for specific condition
    evoked = epochs[condition].average(chs)
    FREQs_temp2 = []
    #This loop iterates for each epoch
    for min, max in FREQ_BANDS:
        FREQs_temp = []
        #psds, freqs = function(epochs[t], fmin=min, fmax=max, bandwidth = 4, picks = EEG_chs)  #PSDs are calculated with this function, giving power values and corresponding frequency bins as output
        psds, freqs = psd_welch(evoked, fmin=min, fmax=max, n_fft = fft_size)
        psds = 10. * np.log10(psds)   #PSDs values are transformed in log scale to compensate for the 1/f natural slope
        for ch in range(len(psds)):
            index_max = np.argmax(np.array(psds[ch][:]))
            freq = freqs[index_max]
            print(index_max) # Should not be zero in all bands (would signify strong 1/f trend)
            FREQs_temp.append(freq)
        if dim_reduc =='avg':
            FREQs_avg = np.average(FREQs_temp)
            FREQs_temp2.append(FREQs_avg)
        if dim_reduc =='mode':
            s = pd.Series(FREQs_temp)
            FREQs_mode = s[s.duplicated()].unique().tolist()
            import itertools

            FREQs_temp2.append(FREQs_mode)
            print(FREQs_temp2)
        if dim_reduc ==None:
            FREQs_temp2.append(FREQs_temp)
    if dim_reduc =='avg':
        FREQs_temp2 = [round(num, 1) for num in FREQs_temp2]
    elif dim_reduc =='mode':
        FREQs_temp2 = list(itertools.chain(*FREQs_temp2))
    FREQs = np.array(FREQs_temp2)
    return FREQs
