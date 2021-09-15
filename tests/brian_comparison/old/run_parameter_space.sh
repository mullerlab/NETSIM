#!/bin/bash

######################
#
# run a parameter search space with both netsim and brian
#
######################

master_param_path=$1
DELAYS=$2

# get current dir
path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# brian 
param_creator_path="/home/tdesbordes/Documents/netrover/param_creator.py"
python_brian_path="$path/brian_net_param_space.py"

# netsim path
binary_path="$(dirname $(dirname $path))/netsim" # get path to netsim binary
# check binary
if ! [[ -f $binary_path ]]; then echo "Incorrect binary path"; exit 1; fi
if ! [[ -x $binary_path ]]; then echo "Binary not excecutable"; exit 1; fi

# NETSIM specific options
P_RELEASE_FLAG=""; EXTERNAL_INPUT_FLAG=""
if [[ "$(basename $binary_path)" == "netsim" ]]; then

	# conditionals for compilation step
	if grep -Fxq "p_release" $1; then
		echo "compiling with p_release..."
		P_RELEASE_FLAG+="RELEASE_PROBABILITY=yes"
	fi

	if grep -Fxq "poisson_rate" $1; then
		echo "compiling with external input..."
		EXTERNAL_INPUT_FLAG+="EXTERNAL_INPUT=yes"
	fi

  # issue MAKE command
	ssh bertha "make -C $(dirname $binary_path) clean; make -C $(dirname $binary_path) $P_RELEASE_FLAG $EXTERNAL_INPUT_FLAG"

fi


# put together path
timestamp=/brian_comp_run_"$(date +%Y%m%d%H%M-%S)"
run_path="$(dirname "$master_param_path")$timestamp"
# generate all parameter files
python $param_creator_path $master_param_path $timestamp

# create filenames
filenames=$(find $run_path -name '*.parameters')
nb_simulations=$(find $run_path -name '*.parameters' | wc -l)

# make output dir
mkdir -p $(dirname $run_path)/outputs

i_simu=0
# simulations loop
for i in ${filenames[@]}; do
  echo "Running simulation nb $i_simu with parameter $i "
  echo "netsim simulation " 
  $binary_path -f $i
  echo "brian simulation " 
  python $python_brian_path $i $i_simu $nb_simulations $DELAYS
  return=$?
  # be sure that the script is ran without error for every parameter file
  while [ $return -ne 0 ]; do
    echo "Error detected when running the python script. Trying to run it again with the same parameters"
    python $python_brian_path $i $i_simu $nb_simulations $DELAYS
    return=$?
  done
#  $binary_file -f $i  | grep 'Rate' >> $output_file
  ((i_simu++))
done

# move brian output
mkdir $(dirname $run_path)/outputs_brian/
mv $run_path/outputs_brian/* $(dirname $run_path)/outputs_brian/
rmdir $run_path/outputs_brian/



# spike analysis binary file
spike_analysis_path="/home/tdesbordes/Documents/ns1_figures/analysis/spike_analysis/spike_analysis"
make -C $(dirname $spike_analysis_path/) clean; make -C $(dirname $spike_analysis_path/)

# Spike analysis
out_path="$(dirname "$master_param_path")"

# remove output files if already exists
rm -f $out_path/fr.txt
rm -f $out_path/cv.txt

ctr=1

for filename in ${filenames[@]}; do
  # check that the output file exists, if not put a nan in fr.txt and cv.txt
  output_path=$(grep -oP '(?<=\$output_path = )[ A-Za-z0-9-/-_]*' $filename)
  output_filecode=$( printf "%08dspk.bin" $ctr )

  if [ ! -f $output_path$output_filecode ]; then
  	echo 'nan' >> $out_path/fr.txt
  	echo 'nan' >> $out_path/cv.txt
  	echo 'not found, adding nan'
  fi
  out=$( ($spike_analysis_path -f $filename) | grep 'average_rate' )
  echo $out | grep -i "average_rate:" | cut -d"," -f2 | cut -d":" -f2 | sed -e 's/^[[:space:]]*//' >> $out_path/fr.txt
  echo $out | grep -i "cv =" | cut -d"," -f4 | cut -d"=" -f2 | sed -e 's/^[[:space:]]*//' >> $out_path/cv.txt

  ctr=$((ctr + 1))
done
     

# run the report creating script
python_report_path="$path/brian_parameter_space_make_report.py"
python $python_report_path $(dirname $run_path)
