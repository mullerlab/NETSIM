#!/bin/bash

python_brian_path=$1
filename=$2
i_simu=$3
nb_simulations=$4
DELAYS=$5

python $python_brian_path $filename $i_simu $nb_simulations $DELAYS
return=$?
# be sure that the script is ran without error for every parameter file
while [ $return -ne 0 ]; do
  echo "Error detected when running the python script. Trying to run it again with the same parameters"
  python $python_brian_path $filename $i_simu $nb_simulations $DELAYS
  return=$?
done
