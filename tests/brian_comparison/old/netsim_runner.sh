#!/bin/bash

######################
#
# rOBSOLETE -- USE run_parameter_space_SGE.sh or run_parameter_space.sh for running both brian and netsim simulations at the same time
#
######################

run_path=$1

path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

binary_path="$(dirname $(dirname $path))/netsim" # get path to netsim binary\

echo $binary_path

# check binary
if ! [[ -f $binary_path ]]; then echo "Incorrect binary path"; exit 1; fi
if ! [[ -x $binary_path ]]; then echo "Binary not excecutable"; exit 1; fi

# get filenames
filenames=$(find $run_path -name '*.parameters')
nb_simulations=$(find $run_path -name '*.parameters' | wc -l)

# simulations loop
for i in ${filenames[@]}; do
  echo "Running simulation with parameter $i"
  $binary_path -f $i
done
