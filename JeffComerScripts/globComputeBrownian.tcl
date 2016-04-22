# jcomer2@uiuc.edu
# Analysis

set analysisScript gridTrajBrownian.tcl
set outDir grid

#set analysisScript currentTraj.tcl
#set outDir current

set simSet [list jehan_80]
set dcdDir [list ../]

set dcdCount 10
set dcdFreq 40000
set startFrame 20
set stride 1
set moiety ds

proc checkExistence {path} {
    if {![file exists $path]} {
	puts "ERROR! The file $path does not exist."
	return 0
    }
    return 1
}

proc lmatch {lst reg} {
    set ret {}
    foreach item $lst {
	if {[regexp $reg $item]} {lappend ret $item}
    }
    return $ret
}

proc lnmatch {lst reg} {
    set ret {}
    foreach item $lst {
	if {![regexp $reg $item]} {lappend ret $item}
    }
    return $ret
}


# Check that the files exist.
puts "\nChecking that the necessary files exist..."
set okay 1
set simDcd {}
set simStruct {}
foreach s $simSet dir $dcdDir {
    # Glob for the dcd list.
    set dcdList [glob -directory $dir "${s}.*.dcd"]
    #set dcdList [lnmatch $dcdList  ".*bp0_.*"]
    set dcdList [lsort $dcdList]
    set dcdList [lrange $dcdList 0 [expr {$dcdCount-1}]]
    
    foreach d $dcdList {
	if {![checkExistence $d]} {set okay 0}
    }
    lappend simDcd $dcdList

    set struct [lindex $dcdList 0]
    if {![checkExistence $struct.pdb]} {set okay 0}
    lappend simStruct $struct
}

if {!$okay} {exit}

# Analyze.
source $analysisScript
foreach s $simSet struct $simStruct dcdList $simDcd {
    puts "Simulation: $s"
    puts "$dcdList"
    compute $s $moiety $struct $dcdList $dcdFreq $outDir $stride
}

exit
