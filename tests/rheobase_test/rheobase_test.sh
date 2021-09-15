#!/bin/bash

###############################
#
# rheobase test for netsim
# $1 = binary file
# $2 = param_creator.py path
#
###############################

# parse inputs
binary_file=$1
param_creator_path=$2

output_file='rheobase_results.txt'
rm -f $output_file # delete the ouput file in case it already exists

# check inputs
if ! [[ -f $binary_file ]]; then echo "Incorrect binary path"; exit 1; fi
if ! [[ -x $binary_file ]]; then echo "Binary not excecutable"; exit 1; fi

# get the directory in which the script is
path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
analytical_nb_spikes_py="$path/analytical_nb_spikes.py"
plot_py="$path/plot_rheobase.py"

# put together path
master_parameter_filename="$path/rheobase_test.master"
timestamp=/rheobase_test_run_"$(date +%Y%m%d%H%M-%S)"
run_path="$(dirname "$master_parameter_filename")$timestamp"
# generate all parameter files
python $param_creator_path $master_parameter_filename $timestamp

# create filenames
filenames=$(find $run_path -name '*.parameters')

# simulations loop
for i in ${filenames[@]}; do
  python $analytical_nb_spikes_py $i >> $output_file
  $binary_file -f $i  | grep 'Rate' >> $output_file
done
python $plot_rheobase $output_file $path
