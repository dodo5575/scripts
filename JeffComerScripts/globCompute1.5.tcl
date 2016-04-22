# jcomer2@uiuc.edu
# Analysis

#set analysisScript basepairLocalDensityTraj1.5.tcl
#set outDir density1.5

set analysisScript basepairLocalEnvTrajBrown.tcl
set outDir env

set moiety dna
set stride 100

set sim {}
#lappend sim {distro_TA.100_run ../dcd/nw_distro_TA.100 nw_distro_TA.100_ ../dcd 5000}
#lappend sim {distro_CG.100_run ../dcd/nw_distro_CG.100 nw_distro_CG.100_ ../dcd 5000}
#lappend sim {distro_TA.1000_run ../dcd/nw_distro_TA.1000 nw_distro_TA.1000_  ../dcd 5000}
#lappend sim {distro_CG.1000_run ../dcd/nw_distro_CG.1000 nw_distro_CG.1000_  ../dcd 5000}
#set want {run1 run2 run3 run4 run5 1cont3 1cont4 1cont5 2cont3 2cont4 2cont5}

lappend sim {distro_TA.100_restrain ../dcd/nw_distro_TA.100 nw_distro_TA.100_  ../dcd 5000}
lappend sim {init_distAggr_TA.100 ../4_brownian/init/init_distAggr_TA.100.brown.0.dcd init_distAggr_TA.100  ../4_brownian/init 5000}
set want {restrain brown.0 brown.1 brown.2 brown.3 brown.4}

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
    foreach {name struct dcdPrefix dcdDir dcdFreq} $s { break }

    #if {![checkExistence $struct.psf]} {set okay 0}
    if {![checkExistence $struct.pdb]} {set okay 0}
    #if {![checkExistence $struct.xsc]} {set okay 0}

    # Glob for the dcd list.
    set dcdList0 [glob -directory $dcdDir "${dcdPrefix}*.dcd"]
    set dcdList {}
    foreach w $want {
	set dcdList [concat $dcdList [lmatch $dcdList0 "$w"]]
    }
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
    foreach {name struct dcdPrefix dcdDir dcdFreq} $s { break }

    puts "Simulation: $s"
    puts "$dcdList"
    compute $name $moiety $struct $dcdList $dcdFreq $outDir $stride
}

exit
