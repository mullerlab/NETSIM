#!/bin/bash

function log {
	echo "`date +"%Y-%m-%d %H:%M:%S"`: $@"
}

function error {
	log "$@" > /dev/stderr
}

# This exists so we can perform error handling in a single place:
function execute_job() {
	$PROGRAM $ARGUMENTS	
	RC=$?
	if [ $RC -ne 0 ]; then
		error "$PROGRAM exited with $RC"
	fi
}

function resume_job() {
	log "Resuming $PROGRAM with $ARGUMENTS"
	execute_job
}

function begin_job() {
	log "Starting $PROGRAM with $ARGUMENTS"
	execute_job
}

# These empty functions exist solely to prevent errors if the user doesn't
# provide their own function:
function pre_run() { true; }
function post_run() { true; }

function run_job() {
	if [ ! -x "$PROGRAM" ]; then
		error "Unable to run: $PROGRAM is not executable!"
		exit 1
	fi

	pre_run
	if [[ -z "$RESTARTED" || "$RESTARTED" -ne "1" ]]; then
		begin_job
	else
		resume_job
	fi
	post_run
}