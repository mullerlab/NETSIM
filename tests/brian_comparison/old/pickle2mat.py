####################
#
# Save all pickle files in the given directory as .mat
#
####################

import os
import sys
from scipy.io import savemat
import pickle


try:
	path = sys.argv[1] 
except:
	path = input("Path not supplied...you can enter it here: ")


pickled_files = sorted( [each for each in os.listdir(path) if each.endswith('.p')] )

for pf in pickled_files:
	f = pickle.load( open(path + pf, 'rb') )
	savemat( path + pf[0:-2] + '.mat', f)

