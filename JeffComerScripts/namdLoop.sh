#! /bin/bash

# Perform a binary search for the right water count.
## called with: script $IBRUN $bindir/namd2 ${cpuList[$i]} ${cpuOffset[$i]} $LOG) &
ibrun=$1
namd=$2
cpus=$3
offset=$3
log=$4
ibrunOpts="VIADEV_RENDEZVOUS_THRESHOLD=50000" ;# ranger specific...

iter=0
maxWater=99652
minWater=$(dc -e "$maxWater 3 * 4 / p")
midWater=$(dc -e "$minWater $maxWater + 2 / p")
targetSize=262
sys=b-dna

#####################################################################
# Evaluate a floating point number expression.
function floatEval()
{
    local stat=0
    local result=0.0
    if [[ $# -gt 0 ]]; then
        result=$(echo "scale=$float_scale; $*" | bc -q 2>/dev/null)
        stat=$?
        if [[ $stat -eq 0  &&  -z "$result" ]]; then stat=1; fi
    fi
    echo $result
    return $stat
}


#####################################################################
# Evaluate a floating point number conditional expression.
function floatCond()
{
    local cond=0
    if [[ $# -gt 0 ]]; then
        cond=$(echo "$*" | bc -q 2>/dev/null)
        if [[ -z "$cond" ]]; then cond=0; fi
        if [[ "$cond" != 0  &&  "$cond" != 1 ]]; then cond=0; fi
    fi
    local stat=$((cond == 0))
    return $stat
}


## to run namd!
runNamd() {
#    $ibrun -n $cpus -o $offset $ibrunOpts $namd $1 &
    #(echo "pretending to run namd"; sleep 1s) &
    
    (echo "runNamd!";) &
    
    (echo $ibrun $namd ++nodelist $TMPDIR/namd-machines +giga +p$NSLOTS $1 > ${1%.namd}.log;) &
    
    $ibrun $namd ++nodelist $TMPDIR/namd-machines +giga +p$NSLOTS $1 > ${1%.namd}.log &
    
    ## keep this here, since it will be used later
    namdJobId=$(jobs -p %) ;# grab the line from jobs corresponding to namd
    echo "namdJobId = $namdJobId"
    #namdJobId=${namdJobId:1:1} pull the number out from it
    sleep .1s
}
killNamd() {
    #kill %$namdJobId
    kill $namdJobId 2> /dev/null
}

createNamdConfig() {
    ## here is where you define what you need to create a new namd file
    ## also created neccessary psf/pdb
    local xstFile=output/run_pore6_${sys}_water${midWater}.xst
    local minTime=10000

    # Check the results of the last simulation.
    if [ -e $xstFile ]; then
	meanSize=$(tclsh meanNptSize.tcl $minTime $xstFile)
	echo "WATERMEANSIZE: $meanSize"

	if floatCond "$meanSize >= 0.0"; then
	    if floatCond "$meanSize < $targetSize"; then
		echo "$meanSize < $targetSize"
		minWater=$midWater
	    else
		echo "$meanSize >= $targetSize"
		maxWater=$midWater
	    fi
	fi
    fi

    # Get the new value of the water count.
    echo "WATERBRACKET: $minWater $maxWater"
    midWater=$(dc -e "$minWater $maxWater + 2 / p")
    echo "WATERCOUNT: $midWater"

    # Check for convergence.
    if (($maxWater - $maxWater)); then
	exit
    fi
    
    # Build the new system.
    vmd -dispdev text -e setSolutionCountRandom.tcl -args $sys $midWater
    vmd -dispdev text -e restrainSiliconWater.tcl -args $sys $midWater

    # Make the new config file.
    ./variableConfig.pl run_pore6_${sys}_water water $midWater
    config=run_pore6_${sys}_water${midWater}.namd
    echo "Created $config"
    ## should probably first check for the existence of output files to obtain the next config file
#    echo $templateConfigFile

    #sed 's/set truncated .*/set truncated $truncated/' $templateConfigFile > $newFile
    #sed 's/set n .*/set n $n/' $templateConfigFile > $newFile
}
jobShouldDie() {
    ## This function runs some diagnostics and 
    
    ## if job should die, return 0

    ## otherwise return -1
    return -1
}
jobIsDead() {
    ## checks to see if namdJobId still exists
    if (( ${#namdJobId} == 0 )); then # namd hasn't run yet, return true
	return 0
    fi

    ## loop through jobs
    jobs > /dev/null ;# this line seems to be needed to update jobs -p
    for pid in $(jobs -p); do
	# echo "$namdJobId == $pid"
	if (( $namdJobId == $pid )); then return -1; fi
    done

    ## no namdJobId not found => return true
    echo jobIsDead!
    return 0
}

truncated=0;
while sleep 1s; do
    ## watch the job
    if jobShouldDie;  then
	killNamd
	createNamdConfig
	runNamd $config

    elif jobIsDead;  then
	createNamdConfig
	runNamd $config
    fi
	
done
