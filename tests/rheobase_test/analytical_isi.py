#from brian2 import *
import sys
import numpy as np
import warnings
warnings.filterwarnings("ignore")

try:
    path = sys.argv[1]
except:
    path = input( 'path to parameter file not provided...you can enter it here: ' )

#path = '/home/tdesbordes/Documents/Salk_internship/netsim/tests/rheobase_test_parameters/01_Ie=0.parameters'

f = open( path, 'r' )
lines = f.readlines()
f.close()

for line in lines:
    try:
        exec(line)
    except:
        continue

# Membrane Parameters
gL = Cm / taum       # = 1e-8 siemens => 100 Mohm
Rm = 1/gL

isi = taum * np.log( (Rm * Ie + El - vreset) / (Rm * Ie + El - vth) ) + taur


print('Ie = ', Ie, '; Theoretical isi (s) = ', isi, end=' ; ')
