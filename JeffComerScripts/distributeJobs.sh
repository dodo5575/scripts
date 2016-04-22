#!/bin/bash
# Author: Jeff Comer <jcomer2@illinois.edu>

# Run jobs on the following machines:
machines=(tbgl-work2 tbgl-work8 tbgl-work11 tbgl-work10 tbgl-work9)
user=jcomer
workDir=/projects/jcomer/protein-dna/cplcRender/frames
runName=./workDistributeRenderTrans.sh

# Get the nubmer of machines.
nMachines=`echo ${machines[*]} | wc -w`
echo "Number of machines: $nMachines"
# Get the number of files.
len=$#
echo "Number of files: $len"
# Split the lists by the number of cores.
let "each = $len/$nMachines + 1"
echo "Files per machine: $each"

mach=0
subList=()
for file in $@; do
   # Add it to the current sublist.
    subList[$count]=$file

   # Increment the count.
    ((count++))

   # We have filled this processor.
    if [ $count -eq $each ]; then
	# Run on this set of files.       
	echo ""
	echo "Running on machine $mach \(${machines[$mach]}\):"
	echo ${subList[*]}

	# Run the command.
	cmd="cd $workDir; $runName ${subList[*]} > log_${mach}.log & disown;"
	sshCmd="ssh -l $user ${machines[$mach]} \"$cmd\""
	#echo $sshCmd
	eval $sshCmd
	
	# Prepare for the next set.
	count=0
	subList=()
	((mach++))
    fi
done

# Run on this set of files.       
echo ""
echo "Running on machine $mach:"
echo ${subList[*]}

# Run the command.
cmd="cd $workDir; $runName ${subList[*]} > log_${mach}.log & disown;"
sshCmd="ssh -l $user ${machines[$mach]} \"$cmd\""
#echo $sshCmd
eval $sshCmd
