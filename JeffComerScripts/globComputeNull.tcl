# jcomer2@uiuc.edu
# Analysis

set analysisScript currentThroughTraj.tcl
set outDir current_null

#set analysisScript cylindricalVelTraj.tcl
#set outDir vel

#set analysisScript cylindricalCurrTraj.tcl
#set outDir curr

#set analysisScript cylindricalCurrentTraj.tcl
#set outDir cyl_current_err

#set analysisScript cylindricalWaterTraj.tcl
#set outDir cyl_water

set sim {}
lappend sim [list "null_pos" pore_no_basepair_nw /scratch/tbgl/jcomer/basepair/dcd 5000]
lappend sim [list "none_pos" ../open/nw_pore_none_basepair ../open/output 5000]
lappend sim [list "no_1M_pos" ../more_ions/nw_pore_no_basepair_1M /projects/jcomer/basepair/more_ions/dcd 5000]
lappend sim [list "no_1M_neg" ../more_ions/nw_pore_no_basepair_1M /projects/jcomer/basepair/more_ions/dcd 5000]

set startFrame 20
set moiety ds

#set simSet [list "null_pos"]
#set simStruct [list pore_no_basepair_nw]
#set dcdDir [list /scratch/tbgl/jcomer/basepair/dcd]

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
foreach s $sim {
    foreach {name struct dir dcdFreq} $s { break }

    if {![checkExistence $struct.psf]} {set okay 0}
    if {![checkExistence $struct.pdb]} {set okay 0}
    if {![checkExistence $struct.xsc]} {set okay 0}

    # Glob for the dcd list.
    set dcdList [glob -directory $dir "nw_*${name}*.dcd"]
    #set dcdList [glob -directory $dir "nw_basepair_no_pos19.dcd"]
    #set dcdList [lnmatch $dcdList  ".*basepairX.*"]
    #set dcdList [lnmatch $dcdList  ".*bp0_.*"]
    set dcdList [lsort $dcdList]
    
    puts "$name [llength $dcdList]"
    if {[llength $dcdList] == 0} { 
	puts "Found no matching files for `$name'!"
	exit 
    }
    
    foreach d $dcdList {
	if {![checkExistence $d]} {set okay 0}
    }
    lappend simDcd $dcdList
}

if {!$okay} {exit}

# Analyze.
source $analysisScript
foreach s $sim dcdList $simDcd {
    foreach {name struct dir dcdFreq} $s { break }

    puts "Simulation: $name"
    #puts "Number of dcd files: [llength $dcdList]"
    compute $name $moiety $struct $dcdList $dcdFreq $outDir $startFrame
}

exit
