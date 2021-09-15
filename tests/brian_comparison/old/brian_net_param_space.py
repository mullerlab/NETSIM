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

# get the parameter file
try:
    path = sys.argv[1]
    i_simu = int(sys.argv[2])
    nb_simulations = int(sys.argv[3])
    DELAYS = sys.argv[4]
except:
    path = input( 'path to parameter file not provided...you can enter it here: ' )

from brian2 import *
from mne import Report
import pickle
import pandas as pd
import os
import time
from scipy.io import savemat

# Read the prameter file
f = open( path, 'r' )
lines = f.readlines()
f.close()

# get the parameters from the parameter file
for line in lines:
    try:
        exec(line)
    except:
        continue

num_saved_param = '' # can be changed to make different versions of the results and average over them
save = True

# make output folder and get the path
if not os.path.isdir( os.path.dirname(path) + os.sep + 'outputs_brian' ):
    os.mkdir( os.path.dirname(path) + os.sep + 'outputs_brian')


save_path = os.path.dirname(path) + os.sep + 'outputs_brian' + os.sep

compute_parameter_space = False
saved_param_string = "netsim_comp_"

we = ge * siemens # excitatory synaptic weight
wi = gi * siemens # inhibitory synaptic weight
del ge, gi

defaultclock.dt = dt * second
t_run = T * second
nb_neurons = N

# Membrane Parameters
C = Cm * farad
taum = taum * second
gL = C / taum       # = 1e-8 siemens => 100 Mohm
EL = El * volt
VT = vth * volt

# Synapse parameters
Ee *= volt
Ei *= volt
taue *= second
taui *= second

# Membrane Equation
eqs = Equations('''
dv/dt = ( gL*(EL-v) + ge*(Ee-v) + gi*(Ei-v) ) * (1./C) : volt (unless refractory)
dge/dt = -ge*(1./taue) : siemens
dgi/dt = -gi*(1./taui) : siemens
''')

print( 'Launching simulation ', i_simu+1, ' / ', nb_simulations, ' ; with ge = ', we / siemens, ' and gi = ', wi / siemens )

# Check whether these parameters haven't already been ran
#if  ( sum((df['we'] == we) & (df['wi'] == wi)) != 0 ) & ( (we != 0) or (wi !=0) ):
#    print( "simulation already done for this parameter set, going to next simulation" )
#    sys.exit()

# Set up the network
P = NeuronGroup(nb_neurons, model=eqs, threshold='v>VT', reset='v=EL', refractory=5*ms, method='euler')
Pe = P[0:int(nb_neurons*8/10)]
Pi = P[int(nb_neurons*8/10)::]

Ce = Synapses( Pe, P, on_pre='ge+=we', multisynaptic_index='Ee' )
Ci = Synapses( Pi, P, on_pre='gi+=wi', multisynaptic_index='Ei' )

Ce.connect( True, p=K/nb_neurons )
Ci.connect( True, p=K/nb_neurons )

## Connection delays

#thetas = linspace( 0, 360, N )
#radius = ( sigma_space * N ) / ( 2 * pi )
#spatial_positions = [ [radius * cos(thetas[i]), radius * sin(thetas[i])] for i in range(N)] # neurons are evenly space on a circle of radius one
#spatial_positions = array(spatial_positions)
#plot(array(spatial_positions)[:,0], array(spatial_positions)[:,1]) # plot the circle

#def arc_length(theta1, theta2):
#    return 2 * pi * radius * (theta2 - theta1) / (360)\
#def edist(a, b): # euclidian distance
#    return sqrt( (b[0] - a[0])**2 + (b[1] - a[1])**2 )
#e_delays = []
#e_conduction_speeds = uniform( min_conduction_speed, max_conduction_speed, len(Pe) )  # meter/s = mm/ms
#for i_syn in range( len(Ce) ) :
#    e_delays.append( arc_length( thetas[Ce.i[i_syn]], thetas[Ce.j[i_syn]] ) / e_conduction_speeds[Ce.i[i_syn]] )
#Ce.delay = array(e_delays) * ms
#del e_delays
#i_delays = []
#i_conduction_speeds = uniform( min_conduction_speed, max_conduction_speed, len(Pi) )  # meter/s = mm/ms
#for i_syn in range( len(Ci) ) :
#    i_delays.append( arc_length( thetas[Ci.i[i_syn]], thetas[Ci.j[i_syn]] ) / i_conduction_speeds[Ci.i[i_syn]] )
#Ci.delay = array(e_delays) * ms
#del i_delays

if DELAYS == 'uniform':
    Ce.delay =  ( uniform( 0, 0.5 / min_conduction_speed, len(Ce) ) + synapse_delay ) * ms
    Ci.delay = ( uniform( 0, 0.5 / min_conduction_speed, len(Ci) ) + synapse_delay ) * ms
elif DELAYS == 'constant':
    Ce.delay = "synapse_delay * second"
    Ci.delay = "synapse_delay * second"



# Initialization = gaussian sampling with mean and sigma given by param file
P.v = np.random.normal( vm_mean, vm_sigma, len(P) ) * volt
P.ge = np.random.normal( ge_mean, ge_sigma, len(P) ) * siemens
P.gi = np.random.normal( gi_mean, gi_sigma, len(P) ) * siemens

# Record the number of spikes
M = SpikeMonitor(P)
Mr = PopulationRateMonitor(P)
#statemon = StateMonitor(P, ('v', 'ge', 'gi'), record=True)


# Run simulation
run( t_run )


def get_ISI_CVs(spikemon):
    ISI_CVs = [] #for every neuron
    spike_trains = spikemon.spike_trains()
    for i_neuron in range(len(spike_trains)):
        ISIs = [] #temporary variable
        for i_spike in range(len(spike_trains[i_neuron])-1):
            ISIs.append(spike_trains[i_neuron][i_spike+1] - spike_trains[i_neuron][i_spike])
        ISI_CVs.append(std(ISIs/ms) / mean(ISIs/ms))
    ISI_CVs = array(ISI_CVs)[logical_not(isnan(ISI_CVs))] # remove NaNs due to totally silent neurons (0/0)
    ISI_CVs = array(ISI_CVs)[logical_not(ISI_CVs==0.0)] # remove zero values due to nearly totally silent neurons (0/some value)
    return ISI_CVs


def check_network_sustained_activity(rates):
    time_limit = int(0.01 / (defaultclock.dt/second/10)) # set a minimum length of time to determine that the netwrok is silent (depends on the simulaton time step)
    rates = rates/Hz
    for i_rate in range(len(rates)):
        try:
            i = 0
            while rates[i_rate+i] == 0.0:# if the network is silent for more than 100time steps
                i += 1
                if i == time_limit:
                    return (i_rate/len(rates)) # it is considered silent starting from this point
        except IndexError: # if not, we return that the activity was sustained all along
            return 1
    return 1


# Save the outputs
output={'we':we/siemens, 'wi':wi/siemens, 'FRs':mean(Mr.rate/Hz), 'sustain_times':check_network_sustained_activity(Mr.rate), 'mean_CVs':mean(get_ISI_CVs(M))}

pickle.dump(output, open(save_path + str(i_simu+1) + '_results_' + saved_param_string + str(num_saved_param) + ".p", 'wb')) ## change i_simu to output_file_code
savemat( save_path + str(i_simu+1) + '_results_' + saved_param_string + str(num_saved_param) + ".mat", output )

# Save all the data 
savemat( save_path + str(i_simu+1) + '_spike_data_' + saved_param_string + str(num_saved_param) + ".mat", {'ids':M.i, 'times':M.t/seconds} )


print('Brian simulation ', i_simu+1, ' completed', '\n')
