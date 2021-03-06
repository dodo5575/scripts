#!/bin/bash
# Author: Jeff Comer <jcomer2@illinois.edu>

#echo "working on: $@"

count=0
for file in $@; do
    if [ $count -eq 0 ]; then
        # Set the output directory.
	outDir=$file
	((count++))
    else
        # Strip the name of the file.
	name=${file##*/}

        # Run meanAutocorrVel.tcl
	./meanAutocorrVel.tcl $outDir/$name $file > $outDir/$name.log
    fi
done
