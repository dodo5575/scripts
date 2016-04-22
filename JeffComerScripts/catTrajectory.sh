#!/bin/bash
# Author: Jeff Comer <jcomer2@illinois.edu>

for first in num*M0.*.dat; do
    prefix=${first%.*}
    name=${prefix%M0*}M
    
# Get the number of this simulation.
    num=${prefix##*.}
    #echo -n "$prefix $num"

    # Open the output the output file.
    outFile=cat_${name}.$num.dat
    echo -n "" > $outFile
    lastTime=0.0

    for f in ${name}*.$num.dat; do
	t=`awk '{print $1}' $f | tail -n1`
	awk "{print \$1+$lastTime,\$2}" $f >> $outFile
	lastTime=$(echo $lastTime + $t | bc -q)
	echo -n "$f "
    done
    echo ""
done
