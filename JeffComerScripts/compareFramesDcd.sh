#!/bin/bash

fileList=`ls -1 ../output/layer*.dcd`

for f in $fileList; do
    name=${f##*/}

    # Check for trajectory with water removed.
    nw=../dcd/nw_${name}
    if [ -e $nw ]; then
	fNum=`catdcd $f | grep "Total frames:" `
	fNum1=${fNum#Total frames:*} 

	nwNum=`catdcd $nw | grep "Total frames:" `
	nwNum1=${nwNum#Total frames:*} 
	
	if [[ $nwNum1 -eq $fNum1 ]]; then
	    echo "# $f OK"
	else
	    echo "$f"
	    echo "# $f and $nw have different numbers of frames!"
	fi
    else 
	echo "$f"
    fi
done
