import mne # Here we import mne, the package that will contain most of the function that we will use today.
from mne.datasets.brainstorm import bst_raw # It is possible to import functions individually. This is helpful since it
                                            # saves time, memory, and makes the calls to the function easier.
from mne.preprocessing import ICA, create_eog_epochs, create_ecg_epochs
from mne.time_frequency import tfr_morlet, psd_multitaper, psd_welch
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat, loadmat
from scipy import stats
import pandas as pd
import os
RT_thresh = 4000

def array_topoplot(toplot, ch_xy, showtitle=False, titles=None, savefig=True, figpath=r'C:\Users\Dell\Jupyter\BrainHackSchool2019_AB\EEG_music_scripts', vmin=-1, vmax=1, mask = None):
    #create fig
    fig, ax = plt.subplots(1,len(toplot), figsize=(10,5))
    #create a topomap for each data array
    for i, data in enumerate(toplot):
        image,_ = mne.viz.plot_topomap(data=data, pos=ch_xy, cmap='Purples', vmin=vmin, vmax=vmax, axes=ax[i], show=False, mask = mask[i, :])
        #option for title
        if showtitle == True:
            ax[i].set_title(titles[i], fontdict={'fontsize': 10, 'fontweight': 'heavy'})
    #add a colorbar at the end of the line (weird trick from https://www.martinos.org/mne/stable/auto_tutorials/stats-sensor-space/plot_stats_spatio_temporal_cluster_sensors.html#sphx-glr-auto-tutorials-stats-sensor-space-plot-stats-spatio-temporal-cluster-sensors-py)
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax[-1])
    ax_colorbar = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(image, cax=ax_colorbar)
    ax_colorbar.tick_params(labelsize=8)
    #save plot if specified
    if savefig == True:
        plt.savefig(figpath, dpi=300)
    plt.show()
    return fig, ax

def compute_t_test(data1, data2):
    results = []
    for freq in range(data1.shape[1]):
        results_temp = []
        for elec in range(data1.shape[2]):
            data1_t_test = data1[:, freq, elec]
            data2_t_test = data2[:,freq, elec]
            results_temp.append(stats.ttest_rel(data1_t_test, data2_t_test))
        results.append(results_temp)
    results = np.array(results)
    t_values = results[:,:,0]
    p_values = results[:,:,1]
    return (results, t_values, p_values)

def p_values_boolean(p_values):
    p_values_boolean = p_values.copy()
    for e in range(p_values.shape[1]):
        for c in range(p_values.shape[0]):
            if p_values[c, e] < 0.05:
                p_values_boolean[c, e] = True
            else:
                p_values_boolean[c, e] = False
    p_values_boolean = np.array(p_values_boolean, dtype='bool')
    return(p_values_boolean)

def taskVSsham(subj, PSD_stage, EPO_stage, FOLDERPATH):
    RUN_LIST = {'pareidolia':['1','2','3','4']}
    task = 'pareidolia'
    task_PSD = []
    sham_PSD = []
    epochs_tot = {}
    for e, run in enumerate(RUN_LIST[task]):
        PSD_name, PSD_path = get_pareidolia_bids(FOLDERPATH, subj, task, run, stage = PSD_stage, cond=None)
        epo_name, epo_path = get_pareidolia_bids(FOLDERPATH, subj, task, run, stage = EPO_stage, cond=None)
        epochs = mne.read_epochs(epo_path)
        epochs_data = epochs.get_data()
        PSD = loadmat(PSD_path)
        PSD = PSD['PSD']
        epochs_tot[run] = epochs
        if e < 3:
            task_PSD.append(PSD)
        else:
            sham_PSD = PSD
    task_PSD = np.array(task_PSD)
    sham_PSD = np.array(sham_PSD)
    return task_PSD, sham_PSD, epochs_tot

def match_n_trials (psds1, psds2):
    if len(psds1)>len(psds2):
        while len(psds1)>len(psds2):
            psds1 = np.delete(psds1, len(psds2), 0)
            #psds1 = np.array(psds1)
            #psds2 = np.array(psds2)
    if len(psds2)>len(psds1):
        while len(psds2)>len(psds1):
            psds2 = np.delete(psds2, len(psds1), 0)
            #psds2 = np.array(psds2)
            #psds1 = np.array(psds1)
    return (psds1, psds2)

def split_FD_anova (data, epochs):
    cond1 = []
    cond2 = []
    cond3 = []
    for e, i in enumerate(epochs.events[:, 2]):
        if i == 70 or i == 770 or i == 7770:
            cond1.append(data[e])
        if i == 71 or i == 771 or i == 7771:
            cond2.append(data[e])
        if i == 72 or i == 772 or i == 7772:
            cond3.append(data[e])
    cond1 = np.array(cond1)
    cond2 = np.array(cond2)
    cond3 = np.array(cond3)
    return cond1, cond2, cond3

def split_by_condition (data, epochs, a, b, condition = None, epochs_type = None):
    print(epochs)
    cond1 = []
    cond2 = []
    if condition == None:
        for e, i in enumerate(epochs.events[:, 2]):
            if i == a:
                cond1.append(data[e])
            if i == b:
                cond2.append(data[e])
        cond1 = np.array(cond1)
        cond2 = np.array(cond2)
    if condition == 'parNOpar':
        cond1 = []
        cond2 = []
        if epochs_type == 'RT':
            for e, i in enumerate(epochs.events[:, 2]):
                if i == 44 or i == 440 or i == 441 or i == 442 or i == 4 or i == 40 or i == 41 or i == 42:
                    cond1.append(data[e])
                if i == 7 or i == 70 or i == 71 or i == 72:
                    cond2.append(data[e])
            cond1 = np.array(cond1)
            cond2 = np.array(cond2)
        if epochs_type == None:
            for e, i in enumerate(epochs.events[:, 2]):
                if i == 7 or i == 70 or i == 71 or i == 72:
                    cond2.append(data[e])
                if i == 77 or i == 770 or i == 771 or i == 772 or i == 777 or i == 7770 or i == 7771 or i == 7772:
                    cond1.append(data[e])
        cond1 = np.array(cond1)
        cond2 = np.array(cond2)
        #cond1, cond2 = match_n_trials(cond1, cond2)
        #data = np.stack((cond1,cond2))
    if condition == 'slowVSfast':
        cond1 = []
        cond2 = []
        if epochs_type == 'RT':
            for e, i in enumerate(epochs.events[:, 2]):
                if i == 44 or i == 440 or i == 441 or i == 442:
                    cond1.append(data[e])
                if i == 4 or i == 40 or i == 41 or i == 42:
                    cond2.append(data[e])
            cond1 = np.array(cond1)
            cond2 = np.array(cond2)
        if epochs_type == None:
            for e, i in enumerate(epochs.events[:, 2]):
                if i == 777 or i == 7770 or i == 7771 or i == 7772:
                    cond1.append(data[e])
                if i == 77 or i == 770 or i == 771 or i == 772:
                    cond2.append(data[e])
            cond1 = np.array(cond1)
            cond2 = np.array(cond2)

    if condition == 'FD01':
        cond1 = []
        cond2 = []
        if epochs_type == 'RT':
            for e, i in enumerate(epochs.events[:, 2]):
                if i == 40 or i == 440:
                    cond1.append(data[e])
                if i == 41 or i == 441:
                    cond2.append(data[e])
            cond1 = np.array(cond1)
            cond2 = np.array(cond2)
        if epochs_type == None:
            for e, i in enumerate(epochs.events[:, 2]):
                if i == 70 or i == 770 or i == 7770:
                    cond1.append(data[e])
                if i == 71 or i == 771 or i == 7771:
                    cond2.append(data[e])
            cond1 = np.array(cond1)
            cond2 = np.array(cond2)

    if condition == 'FD02':
        cond1 = []
        cond2 = []
        if epochs_type == 'RT':
            for e, i in enumerate(epochs.events[:, 2]):
                if i == 40 or i == 440:
                    cond1.append(data[e])
                if i == 42 or i == 442:
                    cond2.append(data[e])
            cond1 = np.array(cond1)
            cond2 = np.array(cond2)
        if epochs_type == None:
            for e, i in enumerate(epochs.events[:, 2]):
                if i == 70 or i == 770 or i == 7770:
                    cond1.append(data[e])
                if i == 72 or i == 772 or i == 7772:
                    cond2.append(data[e])
            cond1 = np.array(cond1)
            cond2 = np.array(cond2)
    if condition == 'FD12':
        cond1 = []
        cond2 = []
        if epochs_type == 'RT':
            for e, i in enumerate(epochs.events[:, 2]):
                if i == 41 or i == 441:
                    cond1.append(data[e])
                if i == 42 or i == 442:
                    cond2.append(data[e])
            cond1 = np.array(cond1)
            cond2 = np.array(cond2)
        if epochs_type == None:
            for e, i in enumerate(epochs.events[:, 2]):
                if i == 71 or i == 771 or i == 7771:
                    cond1.append(data[e])
                if i == 72 or i == 772 or i == 7772:
                    cond2.append(data[e])
            cond1 = np.array(cond1)
            cond2 = np.array(cond2)
    return cond1, cond2


'''def split_by_condition (data, epochs, a, b):
    cond1 = []
    cond2 = []
    for e, i in enumerate(epochs.events[:, 2]):
        if i == a:
            cond1.append(data[e])
        if i == b:
            cond2.append(data[e])
    cond1 = np.array(cond1)
    cond2 = np.array(cond2)
    return cond1, cond2'''

def array_topoplot(toplot, ch_xy, showtitle=False, titles=None, savefig=True, figpath=r'C:\Users\Dell\Jupyter\BrainHackSchool2019_AB\EEG_music_scripts', vmin=-1, vmax=1, mask = None):
    #create fig
    fig, ax = plt.subplots(1,len(toplot), figsize=(10,5))
    #create a topomap for each data array
    for i, data in enumerate(toplot):
        image,_ = mne.viz.plot_topomap(data=data, pos=ch_xy, cmap='Spectral_r', vmin=vmin, vmax=vmax, axes=ax[i], show=False, mask = mask[i, :])
        #option for title
        if showtitle == True:
            ax[i].set_title(titles[i], fontdict={'fontsize': 10, 'fontweight': 'heavy'})
    #add a colorbar at the end of the line (weird trick from https://www.martinos.org/mne/stable/auto_tutorials/stats-sensor-space/plot_stats_spatio_temporal_cluster_sensors.html#sphx-glr-auto-tutorials-stats-sensor-space-plot-stats-spatio-temporal-cluster-sensors-py)
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax[-1])
    ax_colorbar = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(image, cax=ax_colorbar)
    ax_colorbar.tick_params(labelsize=8)
    #save plot if specified
    if savefig == True:
        plt.savefig(figpath, dpi=300)
    plt.show()
    return fig, ax

def compute_t_test(data1, data2):
    results = []
    for freq in range(data1.shape[1]):
        results_temp = []
        for elec in range(data1.shape[2]):
            data1_t_test = data1[:, freq, elec]
            data2_t_test = data2[:,freq, elec]
            results_temp.append(stats.ttest_rel(data1_t_test, data2_t_test))
        results.append(results_temp)
    results = np.array(results)
    t_values = results[:,:,0]
    p_values = results[:,:,1]
    return (results, t_values, p_values)

def get_pareidolia_bids(FOLDERPATH, subj, task , run, stage, cond=None):
    '''
    Constructs BIDS basename and filepath in the SAflow database format.
    '''
    if 'eeg' in stage:
        extension = '.fif'
    if 'preproc' in stage:
        extension = '_raw.fif'
    if 'epo' in stage: # determine extension based on stage
        extension = '.fif'
    if 'PSD' in stage:
        extension = '.mat'
    if 'fig' in stage:
        extension = '.jpg'
    if 'Morlet' in stage:
        extension = '-tfr.h5'
    if 'report' in stage:
        extension = '.html'
    if 'behav' in stage:
        extension = '.csv'
    if 'spectral' in stage:
        extension = '.mat'
    elif 'Complexity' in stage:
        extension = '.mat'
    elif 'sources' in stage:
        extension = '.hd5'
    elif 'events' in stage:
        extension = '.tsv'
    elif 'ARlog' in stage:
        extension = '.hdf5'

    if 'events' in stage:
        pareidolia_bidsname = 'sub-{}_ses-EEG_task-{}_run-{}_{}{}'.format(subj, task, run, stage, extension)
    else:
        if cond == None: # build basename with or without cond
            pareidolia_bidsname = 'sub-{}_ses-EEG_task-{}_run-{}_{}{}'.format(subj, task, run, stage, extension)
        else:
            pareidolia_bidsname = 'sub-{}_ses-EEG_task-{}_run-{}_{}_{}{}'.format(subj, task, run, cond, stage, extension)

    pareidolia_bidspath = os.path.join(FOLDERPATH, 'sub-{}'.format(subj), 'ses-EEG', 'eeg', pareidolia_bidsname)
    return pareidolia_bidsname, pareidolia_bidspath

def reformat_events (events, FDlist, RT_thresh = 4000, task = 'pareidolia', run = None, slowVSfast = False, FD = 'class'):

    if task == 'pareidolia':
        if events[4][2] == 1:
            events[:, 2] = np.where(events[:, 2]==1, 0, events[:, 2])
            events[:, 2] = np.where(events[:, 2]==2, 1, events[:, 2])
            events[:, 2] = np.where(events[:, 2]==3, 2, events[:, 2])
            events[:, 2] = np.where(events[:, 2]==4, 3, events[:, 2])
            events[:, 2] = np.where(events[:, 2]==5, 4, events[:, 2])
        if events[4][2] == 5:
            events[:, 2] = np.where(events[:, 2]==1, 0, events[:, 2])
            events[:, 2] = np.where(events[:, 2]==2, 1, events[:, 2])
            events[:, 2] = np.where(events[:, 2]==3, 2, events[:, 2])
            events[:, 2] = np.where(events[:, 2]==4, 4, events[:, 2])
            events[:, 2] = np.where(events[:, 2]==5, 3, events[:, 2])
        ##Rename Image onset events to inform about pareidolia and FDx`
        i = -1
        RT = []

        if FD == 'all':
            for z in range(len(FDlist)):
                print(FDlist[z])
                FDlist[z] = int(FDlist[z]*10)
                #print(FDlist[z])
        for e in range(len(events)):
            if (events[e][2] == 2 and events[e+1][2] == 4):
                RT.append(events[e+1][0]-events[e][0])

        medianRT = np.median(RT)
        print(medianRT)
        for e in range(len(events)):
            #All image onsets are set to '7'
            if (events[e][2] == 2):
                events[e][2] = 7


            if slowVSfast == True:
                #Late
                if (events[e][2] == 7 and events[e+1][2] == 4 and ((events[e+2][0])-(events[e+1][0]))<=RT_thresh):
                    events[e][2] = 777
                #Early
                if (events[e][2] == 7 and events[e+1][2] == 4 and ((events[e+2][0])-(events[e+1][0]))>RT_thresh):
                    events[e][2] = 77
                # Set FD
                if FD == 'all':
                    if (events[e][2] == 7 or events[e][2] == 77 or events[e][2] == 777):
                        i = i+1
                        events[e][2] = int(str(events[e][2]) + str(FDlist[i]))
                if FD == 'class':
                    FDclass = FD2FDclass(FDlist)
                    #print(len(FDclass))
                    #print(FDclass)
                    if (events[e][2] == 7 or events[e][2] == 77 or events[e][2] == 777):
                        i = i+1
                        print(i)
                        events[e][2] = int(str(events[e][2]) + str(FDclass[i]))
                        #print(events[e][2])


                    if (events[e][2] == 4 and ((events[e+1][0])-(events[e][0]))<=RT_thresh):
                        events[e][2] = int('44'+str(FDclass[i]))
                    if (events[e][2] == 4 and ((events[e+1][0])-(events[e][0]))>RT_thresh):
                        events[e][2] = int('4'+str(FDclass[i]))
            if slowVSfast == False:
                if (events[e][2] == 7 and events[e+1][2] == 4):
                    events[e][2] = 77
                if FD == 'all':
                    if (events[e][2] == 7 or events[e][2] == 77):
                        i = i+1
                        events[e][2] = int(str(events[e][2]) + str(FDlist[i]))
                if FD == 'class':
                    FDclass = FD2FDclass(FDlist)
                    #print(FDclass)
                    if (events[e][2] == 7 or events[e][2] == 77):
                        i = i+1
                        events[e][2] = int(str(events[e][2]) + str(FDclass[i]))
                    if (events[e][2] == 4 and ((events[e+1][0])-(events[e][0]))<=RT_thresh):
                        events[e][2] = int('44'+str(FDclass[i]))
                    if (events[e][2] == 4 and ((events[e+1][0])-(events[e][0]))>RT_thresh):
                        events[e][2] = int('4'+str(FDclass[i]))
        '''for e in range(len(events)):
            if (events[e][2] == 7 or events[e][2] == 70 or events[e][2] == 71 or events[e][2] == 72):# and events[e+1][2] != 4):
                print('found_new')
                new_timing = events[e][0] + medianRT
                events = np.insert(events, e+1, np.array((new_timing, 0, 55)), 0)'''
    if task == 'RS':
        if run == '1':
            events[:, 2] = np.where(events[:, 2]==1, 100, events[:, 2])
            events[:, 2] = np.where(events[:, 2]==2, 101, events[:, 2])
        if run == '2':
            events[:, 2] = np.where(events[:, 2]==1, 200, events[:, 2])
            events[:, 2] = np.where(events[:, 2]==2, 201, events[:, 2])

    return events, medianRT

def FD2FDclass (FDlist):
    FDclass = FDlist.copy()
    for i in range(len(FDlist)):
        if FDlist[i] <= 1.1:
            FDclass[i] = 0
        if FDlist[i] >= 1.6:
            FDclass[i] = 2
        if FDlist[i] > 1.1 and FDlist[i] < 1.6:
            FDclass[i] = 1
    return (FDclass)

def arrange_dataframe(all_frame):
  #Create n_objets without zeros
    all_frame['positive_n_objets'] = all_frame['n_objets.response']
  # Creates a boolean column that represents the presence or absence of pareidolia
    all_frame['boolean_n_objets'] = all_frame['n_objets.response']
    all_frame.loc[all_frame['boolean_n_objets']== 'None',['boolean_n_objets']] = 0
    all_frame.loc[all_frame['boolean_n_objets']== '2',['boolean_n_objets']] = '1'
    all_frame.loc[all_frame['boolean_n_objets']== '3',['boolean_n_objets']] = '1'
    all_frame.loc[all_frame['boolean_n_objets']== '4',['boolean_n_objets']] = '1'
    all_frame.loc[all_frame['boolean_n_objets']== '5',['boolean_n_objets']] = '1'
    # Transforms 'none' in '0'
    all_frame.loc[all_frame['n_objets.response']== 'None',['n_objets.response']] = 0
    # Creates a 'bloc' column
    all_frame.loc[all_frame['trials.thisRepN']== 0,['trials.thisRepN']] = 1
    all_frame.loc[all_frame['trials_2.thisRepN']== 0,['trials.thisRepN']] = 2
    all_frame.loc[all_frame['trials_3.thisRepN']== 0,['trials.thisRepN']] = 3
    all_frame.loc[all_frame['trials_sham.thisRepN']== 0,['trials.thisRepN']] = 4
    all_frame['bloc'] = all_frame['trials.thisRepN']
    # Delete rows from practice bloc
    all_frame.rename(columns={'trials_4.thisRepN':'practice'}, inplace=True)
    all_frame = all_frame.drop(all_frame[all_frame.practice == 0].index)
    # Change 'objects' to 'floats'
    all_frame['boolean_n_objets'] = all_frame['boolean_n_objets'].astype('float')
    all_frame['n_objets.response'] = all_frame['n_objets.response'].astype('float')
    all_frame.loc[:, 'vividness.response'] = all_frame.loc[:, 'vividness.response'].replace('None', np.NaN)
    all_frame['vividness.response'] = all_frame['vividness.response'].astype('float')
    all_frame.loc[:, 'vividness.rt'] = all_frame.loc[:, 'vividness.rt'].replace('None', np.NaN)
    all_frame['vividness.rt'] = all_frame['vividness.rt'].astype('float')
    all_frame.loc[:, 'positive_n_objets'] = all_frame.loc[:, 'positive_n_objets'].replace('None', np.NaN)
    all_frame['positive_n_objets'] = all_frame['positive_n_objets'].astype('float')

    all_frame.loc[:, 'reaction_time'] = all_frame.loc[: , 'reaction_time.rt']
    all_frame.loc[:, 'reaction_time'] = all_frame.loc[:, 'reaction_time'].clip(0, 8)

    return all_frame
