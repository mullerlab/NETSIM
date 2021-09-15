#!/bin/bash

###############################
#
# rheobase test for netsim -- need to run netsim before
# $1 = param_creator.py path
# $2 = binary spike analysis file
#
###############################

run_path=$1

# parse inputs
binary_spike_analysis_file="/home/tdesbordes/Documents/ns1_figures/analysis/spike_analysis/spike_analysis"

output_file='rheobase_isi_results.txt'
rm -f $output_file # delete the ouput file in case it already exists


# get the directory in which the script is
path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
analytical_isi_py="$path/analytical_isi.py"

# create filenames
filenames=$(find $run_path -name '*.parameters')

# simulations loop
for i in ${filenames[@]}; do
  python $analytical_isi_py $i >> $output_file
  $binary_spike_analysis_file -f $i  | grep -o 'mean_isi = ........' >> $output_file
done
