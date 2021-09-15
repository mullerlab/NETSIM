#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 15:14:26 2018

@author: tdesbordes
"""
import matplotlib as mpl
mpl.use('Agg') # to able to be run on remote machine
import sys
from numpy import array
txt_path = sys.argv[1]
save_path = sys.argv[2]

txt_file = open( txt_path )
output = txt_file.readlines()
txt_file.close()

Ie = []
th_nb_spikes = []
actual_nb_spikes = []

for line in output:
    Ie.append( float( line.split(';')[0].split('=')[1] ))
    th_nb_spikes.append( float( line.split(';')[1].split('=')[1] ))
    actual_nb_spikes.append( float( line.split(';')[2].split(':')[1] ))

#for i, el in enumerate(output.split(';')):
#    if i % 3 == 0:
#        Ie.append( float(el.split('=')[1]) )
#    elif i % 3 == 1:
#        th_nb_spikes.append( float(el.split('=')[1]) )
#    else:
#        actual_nb_spikes.append( float(el.split(':')[1]) )
    
import matplotlib.pyplot as plt

fig = plt.figure()
plt.plot(Ie, th_nb_spikes, c='r', markersize=3, label='Theoretical')
plt.plot(Ie, actual_nb_spikes, '--', c='b', markersize=3, label='Simulation')
plt.xlabel('Input current (A)')
plt.ylabel('Firing rate (Hz)')
plt.legend()
fig.savefig(save_path+'rheobase_plot.png')

print('Saved figure comparing theoretical and simulation firing rate')

difference = abs( sum(th_nb_spikes) / sum(actual_nb_spikes) - 1 )
print('The sum of all firing rates differ by about ', difference, '%')