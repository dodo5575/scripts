# jcomer2@uiuc.edu
# Analysis

set analysisScript currentThroughTraj.tcl
set outDir current

set moiety dna

set sim {}
lappend sim [list "triplet_cgg" ../nw_triplet_cgg_ions ../nw 5000 200]
lappend sim [list "triplet_gcc" ../nw_triplet_gcc_ions ../nw 5000 200]

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
    foreach {name struct dir dcdFreq startFrame} $s { break }

    if {![checkExistence $struct.psf]} {set okay 0}
    if {![checkExistence $struct.pdb]} {set okay 0}
    if {![checkExistence $struct.xsc]} {set okay 0}

    # Glob for the dcd list.
    set dcdList [glob -directory $dir "nw_*${name}*.dcd"]
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
    foreach {name struct dir dcdFreq startFrame} $s { break }

    puts "Simulation: $name"
    puts "Number of dcd files: [llength $dcdList]"
    compute $name $moiety $struct $dcdList $dcdFreq $outDir $startFrame
}

exit
