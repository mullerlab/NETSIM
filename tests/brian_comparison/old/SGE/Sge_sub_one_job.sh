#!/bin/bash
#$ -S /bin/bash
#$ -cwd

source SGE/Sge-common.sh

binary_file=$1
param_file=$2


# This should be a fully-qualified path to the executable:
PROGRAM=$binary_file
ARGUMENTS="-f $param_file"

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
		function pre_run {
      # get output path
      output_path=$(grep -oP '(?<=\$output_path = )[ A-Za-z0-9-/-_]*' $param_file)
      # Create temp output folder
      new_folder=$(mktemp -d --dry-run)
      new_folder=$(basename $new_folder)
      new_folder=$(grep -oP '(?<=tmp.)[ A-Za-z0-9-/]*' <<< $new_folder)
      mkdir "/scratch/$new_folder"

      # change the output path to write in the scratch directory
      sed -i "\|\$output_path | c \$output_path = /scratch/$new_folder/" $param_file
		}

		function post_run {
      # copy files to output directory using rsync
	rsync -a "/scratch/$new_folder/" $output_path
      # delete tmp dir
      rm -rf "/scratch/$new_folder/"
      # restaure correct $output_path in the parameter file
      sed -i "\|\$output_path | c \$output_path = $output_path" $param_file

		}

run_job
