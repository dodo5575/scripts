# jcomer2@uiuc.edu
# Analysis

set analysisScript basepairLocalDensityTraj.tcl
set outDir density

set dcdDir ../dcd/
set moiety dna
set stride 1

set sim {}
# lappend sim {distro_TA.100_run ../dcd/nw_distro_TA.100 distro_TA.100_run 5000}
# lappend sim {distro_CG.100_run ../dcd/nw_distro_CG.100 distro_CG.100_run 5000}
#lappend sim {distAggr_TA.100_run ../dcd/nw_distro_TA.100 distAggr_TA.100_run 9984}
#lappend sim {distAggr_CG.100_run ../dcd/nw_distro_CG.100 distAggr_CG.100_run 9984}

#lappend sim {distro_TA.1000_run ../dcd/nw_distro_TA.1000 distro_TA.1000_run 5000}
#lappend sim {distro_CG.1000_run ../dcd/nw_distro_CG.1000 distro_CG.1000_run 5000}
lappend sim {distAggr_TA.1000_run ../dcd/nw_distro_TA.1000 distAggr_TA.1000_run 9984}
lappend sim {distAggr_CG.1000_run ../dcd/nw_distro_CG.1000 distAggr_CG.1000_run 9984}

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
    foreach {name struct dcdPrefix dcdFreq} $s { break }

    if {![checkExistence $struct.psf]} {set okay 0}
    if {![checkExistence $struct.pdb]} {set okay 0}
    #if {![checkExistence $struct.xsc]} {set okay 0}

    # Glob for the dcd list.
    set dcdList [glob -directory $dcdDir "nw_*${name}*.dcd"]
    set dcdList [lnmatch $dcdList "0\.dcd"]
    set dcdList [lsort $dcdList]
    
    foreach d $dcdList {
	if {![checkExistence $d]} {set okay 0}
    }
    lappend simDcd $dcdList
}

if {!$okay} {exit}

# Analyze.
source $analysisScript
foreach s $sim dcdList $simDcd {
    foreach {name struct dcdPrefix dcdFreq} $s { break }

    puts "Simulation: $s"
    puts "$dcdList"
    compute $name $moiety $struct $dcdList $dcdFreq $outDir $stride
}

exit
