fileList=`ls -1S | head -n15`

for f in $fileList; do
    # Check for trajectory with water removed.
    nw=nw_${f}
    if [ -e $nw ]; then
	fNum=`catdcd $f | grep "Total frames:" `
	fNum1=${fNum#Total frames:*} 

	nwNum=`catdcd $nw | grep "Total frames:" `
	nwNum1=${nwNum#Total frames:*} 
	
	if [[ $nwNum1 -eq $fNum1 ]]; then
	    prefix=${f%.*}
	    outFile=$prefix.1ns.dcd
	    rm $outFile
	    catdcd -o $outFile -stride 20 $f
	    rm $f
	else
	    echo "$f and $nw have different numbers of frames!"
	fi
    fi
done