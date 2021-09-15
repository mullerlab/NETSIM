#!/bin/bash

######################
#
# run a parameter search space with both netsim and brian
#
######################

master_param_path=$1
DELAYS=$2

param_creator_path="~/Documents/netrover/param_creator.py"

path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# netsim path
binary_path="$(dirname $(dirname $path))/netsim" # get path to netsim binary
sge_path_netsim="$path/SGE/Sge_sub_one_job.sh"
# check binary
if ! [[ -f $binary_path ]]; then echo "Incorrect binary path"; exit 1; fi
if ! [[ -x $binary_path ]]; then echo "Binary not excecutable"; exit 1; fi

# NETSIM specific options
P_RELEASE_FLAG=""; EXTERNAL_INPUT_FLAG=""
if [[ "$(basename $binary_file)" == "netsim" ]]; then

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
	ssh bertha "make -C $(dirname $2) clean; make -C $(dirname $2) $P_RELEASE_FLAG $EXTERNAL_INPUT_FLAG"

fi

# brian paths
# path to script to submit job to the grid engine
sge_path_brian="$path/SGE/Sge_sub_one_job_brian.sh"
run_job_path_brian="$path/SGE/Sge_run_one_job_brian.sh"

python_brian_path="$path/brian_net_param_space.py"

# put together path
# master_param_path="$path/brian_comparison_test.master"
timestamp=/brian_comp_run_"$(date +%Y%m%d%H%M-%S)"
run_path="$(dirname "$master_param_path")$timestamp"
# generate all parameter files
ssh bertha "python $param_creator_path $master_param_path $timestamp"

# create filenames
filenames=$(find $run_path -name '*.parameters')
nb_simulations=$(find $run_path -name '*.parameters' | wc -l)

# make output dir
mkdir -p $(dirname $run_path)/outputs

i_simu=0
# simulations loop
for filename in ${filenames[@]}; do
  qsub -l mf=130G $sge_path_brian $run_job_path_brian $python_brian_path $filename $i_simu $nb_simulations $DELAYS
  qsub -l mf=10G $sge_path_netsim $binary_file $filename
  (( i_simu++ ))
done

# run the report creating script
#python_report_path="$path/brian_parameter_space_report.py"
#python $python_report_path $run_path
