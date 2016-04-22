#/bin/bash
# Author: Jeff Comer <jcomer2@illinois.edu>

fileList=instance_short/*.dat
outDir=correlation_short
cores=7

# Get the number of files.
len=`echo $fileList | wc -w`
echo "Number of files: $len"
# Split the lists by the number of cores.
let "each = $len/$cores + 1"
echo "Files per core: $each"

proc=0
count=0
subList=()
for file in $fileList; do
    # Add it to the current sublist.
    subList[$count]=$file

     # Increment the count.
    ((count++))
    
    # We have filled this processor.
    if [ $count -eq $each ]; then
	# Run on this set of files.       
	echo ""
	echo "Running on core $proc:"
	echo ${subList[*]}
	./workMeanAutocorrVel.sh  $outDir ${subList[*]} & disown

	# Prepare for the next set.
	count=0
	subList=()
	((proc++))
    fi
done

if [ $count -ne 0 ]; then
  # Run on this set of files.       
	echo ""
	echo "Running on core $proc:"
	echo ${subList[*]}
	./workMeanAutocorrVel.sh $outDir ${subList[*]} & disown
fi
