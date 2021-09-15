  #!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e error-sge
#$ -o output-sge

source SGE/Sge-common.sh

run_job_path=$1
python_brian_path=$2
filename=$3
i_simu=$4
nb_simulations=$5
DELAYS=$6


# This should be a fully-qualified path to the executable:
PROGRAM=$run_job_path
ARGUMENTS="$python_brian_path $filename $i_simu $nb_simulations $DELAYS"

#"$param_path $param_file $num_param_file"

# OPTIONAL: Define other variables needed by your program:
# PATH=$HOME/bin:$HOME/projects/my_program
# DISPLAY=myworkstation:0.0 # Display graphs on my workstation

# OPTIONAL: Defined a function which customizes the behaviour of your
#	job by running additional commands at certain stages:
# 	begin_job()		Called to start your program normally
# 	resume_job() 	Called to restart your program after it has been
# 								checkpointed or resumed
# 	pre_run()			Called before a job is started or resumed
# 	post_run()		Called after a job has exited
#
# EXAMPLES
#
# 	function resume_job {
# 		$PROGRAM --resume=/path/to/my/checkpoint_file $ARGUMENTS
# 	}
#
#		function pre_run {
#			rsync -a $HOME/projects/big_data_1/ $TMPDIR/big_data_1/
#		}
#		function post_run {
#			rsync -a $TMPDIR/big_data_1/ $HOME/projects/big_data_1/
#		}

run_job
