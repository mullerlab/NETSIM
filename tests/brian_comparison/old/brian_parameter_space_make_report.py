#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 14:45:48 2018

@author: tdesbordes
"""
import matplotlib as mpl # to enable to run on the cluster
mpl.use('Agg') # to able to be run on remote machine
import warnings
warnings.filterwarnings("ignore")
import sys

# get the path to the parameter folder to save the report
try:
    path = sys.argv[1]
except:
    path = input( 'path to folder not provided...you can enter it here: ' )

from brian2 import *
from mne import Report
import pickle
import pandas as pd
from scipy.io import loadmat
import os
import re

N = 10000

nb_runs = 1

report_name = "Sustained_activity_report_netsimcomp_ave.html"
figsize = 1.5 # make big figures
saved_param_string = "results_netsim_comp_"    

## Helper functions
def get_ISI_CVs(mdat):
    ids = mdat['ids']
    ids = [ ids[i][0] for i in range(len(ids)) ]
    times = mdat['times']
    times = [ times[i][0] for i in range(len(times)) ]
    nb_neurons = int( max( ids ) - 1 )
    ISI_CVs = [] #for every neuron
    # spike_trains = mdat['ids']
    for i_neuron in range( nb_neurons ):
        ISIs = [] #temporary variable
        spike_times_for_this_neuron = array(times)[array(ids) == i_neuron + 1]
        for i_spike in range( len(spike_times_for_this_neuron) - 1 ):
        # for i_spike in range(len(spike_trains[i_neuron])-1):
            ISIs.append(spike_times_for_this_neuron[i_spike+1] - spike_times_for_this_neuron[i_spike])
            # ISIs.append(spike_trains[i_neuron][i_spike+1] - spike_trains[i_neuron][i_spike])
        ISI_CVs.append(std(ISIs/ms) / mean(ISIs/ms))
    ISI_CVs = array(ISI_CVs)[logical_not(isnan(ISI_CVs))] # remove NaNs due to totally silent neurons (0/0)
    ISI_CVs = array(ISI_CVs)[logical_not(ISI_CVs==0.0)] # remove zero values due to nearly totally silent neurons (0/some value)
    return ISI_CVs

def save_fig_to_report(fig, caption, section, figsize=figsize):
    global report
    report.add_figs_to_section(fig, captions=caption, section=section, scale=figsize)


def plot_2D_parameter_space(df, fr_vmin=0, fr_vmax=50, cv_vmin=0, cv_vmax=2):
    import seaborn as sns
    sns.set(font_scale=2)  # crazy big
    # FRs
    fig1, ax1 = subplots(figsize=(9, 7))
    pivot = df.pivot(index='we', columns='wi', values='FRs')
    sns.heatmap(pivot, ax=ax1, vmin=fr_vmin, vmax=fr_vmax)
    ax1.invert_yaxis()
    ax1.collections[0].colorbar.set_label('Firing rate (Hz)')
    xticks = list(ax1.get_xticklabels())
    yticks = list(ax1.get_yticklabels())
#    coord_rect = [WIs.index(wi), WEs.index(we)]
    for ix in range(len(xticks)):
        if ix  % 4 != 0:
            xticks[ix] = ''
    ax1.set_xticklabels(xticks)
    for iy in range(len(yticks)):
        if iy  % 4 != 0:
            yticks[iy] = ''
    ax1.set_yticklabels(yticks)
#    ax1.add_patch(Rectangle(coord_rect, 1, 1, fill=False, edgecolor='blue', lw=2))
    plt.tight_layout()
    title('Firing rate parameter search')#'Distribution of firing rates depending on the strength of inhibition and excitation')
#    # Sustain times
#    fig2, ax2 = subplots(figsize=(9, 7))
#    pivot = df.pivot(index='we', columns='wi', values='sustain_times')
#    sns.heatmap(pivot, ax=ax2)
#    ax2.invert_yaxis()
#    ax2.collections[0].colorbar.set_label('Time of sustain activity (s)')
#    ax2.set_xticklabels(xticks)
#    ax2.set_yticklabels(yticks)
##    ax2.add_patch(Rectangle(coord_rect, 1, 1, fill=False, edgecolor='blue', lw=2))
#    plt.tight_layout()
#    title('Sustaining time parameter search')#'Distribution of the network sustaining time depending on the strength of inhibition and excitation')
    # mean_CVs
    fig3, ax3 = subplots(figsize=(9, 7))
    pivot = df.pivot(index='we', columns='wi', values='mean_CVs')
    sns.heatmap(pivot, ax=ax3, vmin=cv_vmin, vmax=cv_vmax)
    ax3.invert_yaxis()
    ax3.collections[0].colorbar.set_label('Coefficient of variation')
    ax3.set_xticklabels(xticks)
    ax3.set_yticklabels(yticks)
#    ax3.add_patch(Rectangle(coord_rect, 1, 1, fill=False, edgecolor='blue', lw=2))
    plt.tight_layout()
    title('Coefficient of variation parameter search')#'Distribution of the mean variation coefficient depending on the strength of inhibition and excitation')
    return fig1, fig3

def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]


# for i_run in arange(nb_runs) + 1:
i_run = 1
    
brian_output_folder = path + '/outputs_brian/' #+ str(i_run) + '/'
netsim_output_folder = path + '/outputs/' #+ str(i_run) + '/'
# netsim_fr_path = path + '/fr.txt' #+ str(i_run) + '.mat'
# netsim_cv_path = path + '/cv.txt'
    
# import and put together brian data
pickled_files = sorted( [each for each in os.listdir(brian_output_folder) if each.endswith('.p')] )
# sort the pickle files
pickled_files.sort(key=natural_keys)
#init dict
grouped = {'we':[], 'wi':[], 'FRs':[], 'mean_CVs':[], 'sustain_times':[]}
# loop to fill dict
for pick in pickled_files:
    data = pickle.load( open(brian_output_folder + pick, "rb") )
    for key in grouped.keys():
        grouped[key].append( data[key] )
     
# if first run, then create the array, if not, just add the values
if i_run == 1:
    df_brian = pd.DataFrame.from_dict(grouped)
else:
    df_brian = df_brian + pd.DataFrame.from_dict(grouped)

# import netsim data
mat_files = sorted( [each for each in os.listdir(netsim_output_folder) if each.endswith('.mat')] )# file = open( netsim_fr_path, "r")
# netsim_fr = file.readlines()
# file.close()
# file = open( netsim_cv_path, "r")
# netsim_cv = file.readlines()
# file.close()
netsim_fr, netsim_cv = [], []
for mfile in mat_files:
    mdat = loadmat( netsim_output_folder + mfile )
    netsim_fr.append( len(mdat['times']) / N )
    print( netsim_fr[-1] )
    netsim_cv.append( mean( get_ISI_CVs( mdat ) ) )
    print( netsim_cv[-1] )


# # change string to numeric
# for i in range( len(netsim_fr) ):
#     netsim_fr[i] = float( netsim_fr[i][0:-1] )
#     netsim_cv[i] = float( netsim_cv[i][0:-1] )


# if first run, then create the array, if not, just add the values
if i_run == 1:
    # transform the arrays into pandas dataframe
    df_netsim = df_brian.copy()
    df_netsim['FRs'] = netsim_fr
    df_netsim['mean_CVs'] = netsim_cv
else:
    df_netsim['FRs'] += netsim_fr
    df_netsim['mean_CVs'] += netsim_cv

# divide the values by the number of runs to make the average over runs
df_brian = df_brian / nb_runs
df_netsim = df_netsim / nb_runs

# make figs for brian
fig1_brian, fig3_brian = plot_2D_parameter_space(df_brian)
# make figs for netsim
fig1_netsim, fig3_netsim = plot_2D_parameter_space(df_netsim)

# subtract the 2 dataframes
df_diff = df_brian - df_netsim
# restore the we and wi in the new dataframe
df_diff['we'] = df_brian['we'] ; df_diff['wi'] = df_brian['wi']

# make figs for the difference between the two
fig1_diff, fig3_diff = plot_2D_parameter_space(df_diff, fr_vmin=-10, fr_vmax=10, cv_vmin=-1, cv_vmax=1)


report = Report()

save_fig_to_report([fig1_brian], caption=['Brian\' Firing rates'], section='FRs', figsize=figsize)
save_fig_to_report([fig1_netsim], caption=['Netsim\' Firing rates'], section='FRs', figsize=figsize)
save_fig_to_report([fig1_diff], caption=['Brian - Netsim\' Firing rates'], section='FRs', figsize=figsize)

#save_fig_to_report([fig2_brian], caption=['Brian\'s sustaining time'], section='STs', figsize=figsize)
#save_fig_to_report([fig2_netsim], caption=['Netsim\'s sustaining time'], section='STs', figsize=figsize)
#save_fig_to_report([fig2_diff], caption=['Brian - Netsim\'s sustaining time'], section='STs', figsize=figsize)

save_fig_to_report([fig3_brian], caption=['Brian\'s coefficient of variation'], section='CVs', figsize=figsize)
save_fig_to_report([fig3_netsim], caption=['Netsim\'s coefficient of variation'], section='CVs', figsize=figsize)
save_fig_to_report([fig3_diff], caption=['Brian - Netsim\'s coefficient of variation'], section='CVs', figsize=figsize)


# Save the report
report.save(path + report_name, overwrite=True, open_browser=False)
