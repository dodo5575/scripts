#!/bin/bash
# Author: Jeff Comer <jcomer2@illinois.edu>

lastRun=0
nextRun=1

for config in cubo*.brown; do
    prefix=${config%.*}
    
    last=${prefix}${lastRun}
    next=${prefix}${nextRun}

	# Run through all restart files with this name.
    for res in $last.*.restart; do
	    # Does the last restart file exist?
	if [ ! -f $res ]; then
	    echo "Warning: The pattern $res does not exist... Skipping"
	    continue
	fi
	
	    # Get the number of this simulation.
	pre=${res%.*}
	num=${pre##*.}

	    # Run the next simulation if there is the last restart file exists,
	    # but no log file for the next simulation.
	nextLog=$next.$num.log
	if [ ! -f $nextLog ]; then
	    echo $next.$num
	    ./condorSubmission1.pl $config $next.$num $last.$num 1
	    condor_submit $next.$num.sub
	else
	    echo "Warning: $nextLog already exists... Skipping"
	fi
    done
    
done
