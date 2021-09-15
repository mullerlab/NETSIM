#!/home/tdesbordes/anaconda3/envs/general3/bin/python

#  Brian2 networks for: Signal Propagation and Logic Gating in Networks of Integrate-and-Fire Neurons, Vogels and Abbott 2005
# Théo Desbordes
# 26/09/2017
# USe on python 3.6 with cythons

from brian2 import *
from mne import Report
import pickle
import pandas as pd
import sys
import os

#try:
#    path = sys.argv[1]
#except:
#    path = input( 'path to parameter file not provided...you can enter it here: ' )

#path = "/home/tdesbordes/Documents/Salk_internship/netsim/tests/rheobase_test/rheobase_test_run_201803201624-42/03_Ie=2e-10.parameters"

folder_path = "/home/tdesbordes/Documents/Salk_internship/netsim/tests/rheobase_test/rheobase_test_run_201803201624-42/"

## lop over all files 
for path in os.listdir( folder_path )[2:4]:
    
    if path == 'lookup_table.csv':
        continue
    
    f = open( folder_path + path, 'r' )
    lines = f.readlines()
    f.close()
    
    for line in lines:
        try:
            exec(line)
        except:
            continue
    
    
    ## theoretical rates 
    # Membrane Parameters
    gL = Cm / taum       # = 1e-8 siemens => 100 Mohm
    Rm = 1/gL
    
    isi = taum * np.log( (Rm * Ie + El - vreset) / (Rm * Ie + El - vth) ) + taur
    rate = 1 / isi
    
    if np.isnan(rate):
        rate = 0
    
    
    defaultclock.dt = dt * second
    t_run = T * second
    nb_neurons = N
    
    # Membrane Parameters
    C = Cm * farad
    taum = taum * second
    gL = C / taum    
    Rm = 1/gL
    EL = El * volt
    VT = vth * volt
    
    # Synapse parameters
    Ee *= volt
    Ei *= volt
    taue *= second
    taui *= second
    Ie *= amp
    vreset *= volt
    
    # Membrane Equation
    eqs = Equations('''
    dv/dt = ( (EL-v) + Ie * (1./gL) ) / taum : volt (unless refractory)
    ''')
    
    # Set up the network
    P = NeuronGroup(nb_neurons, model=eqs, threshold='v>VT', reset='v=vreset', refractory=5*ms, method='exact')
    
    we = 3.327602336e-09 * siemens # excitatory synaptic weight
    wi = 8e-08 * siemens # inhibitory synaptic weight
    
    #Ce = Synapses(Pe, P, on_pre='ge+=we' )
    #Ce.connect( True, p=200/nb_neurons )
    #Ci = Synapses(Pi, P, on_pre='gi+=wi' )
    #Ci.connect( True, p=200/nb_neurons )
    
    P.v = vr * volt
    
    
    # Record the number of spikes
    M = SpikeMonitor(P)
    Mr = PopulationRateMonitor(P)
    statemon = StateMonitor(P, ('v'), record=True)
    
    # Store the initial state
    store()
    
    # Run simulation
    run( t_run )
    
    
    
    # print result
    print('Ie = ', Ie, '; Theoretical firing rate (Hz) = ', rate, '  ;  Simulation rate (Hz) = ', mean(Mr.rate))



## Plots
def plot_raster(spikemon, nb_neuron='all', save=True):
    fig = figure()
    if nb_neuron == 'all':
        plot(spikemon.t/second, spikemon.i, '.k', markersize=0.8)
    else:
        indexes = linspace(0,nb_neuron-1,nb_neuron)
        i = []
        t = []
        for index in indexes:
            i = concatenate((i, spikemon.i[spikemon.i == index]))
            t = concatenate((t, spikemon.t[spikemon.i == index]))
        plot(t, i, '.k')
    xlabel('Time (s)')
    ylabel('Neurons')
    frame = plt.gca()
    frame.axes.set_yticks([])
    # frame.axes().get_yaxis().set_ticks([])
    return fig

def plot_average_FR(spikemon):
    FR = []
    for i in range(len(P)):
        FR.append(sum(spikemon.i == i))
    bins = linspace(0,100,100)
    fig = figure()
    hist(FR, bins=bins)
    vlines(mean(FR), 0, 1000)
    xlabel("Average firing rate (Hz)")
    ylabel("Number of neurons")
    return fig

def plot_overall_FR(ratemon):
    FR = (ratemon.rate)/Hz
    fig = figure()
    plot(ratemon.t/second, FR)
    xlabel("Time (s)")
    ylabel("Overall firing rate (Hz)")
    ylim(-10, max(FR)+10)
    return fig


fig1 = plot_raster(M)
fig2 = plot_raster(M, 250)
#fig3 = plot_average_FR(M)
#fig4= plot_overall_FR(Mr)
