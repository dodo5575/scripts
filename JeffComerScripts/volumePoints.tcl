# Author: Jeff Comer <jcomer2@illinois.edu>

set points 60
set selText "segname ADNA"
# Input:
set psf just_cyt_pot.psf
set pdb just_cyt_pot.pdb
# Output:
set outFile surface.dat

source $env(HOME)/scripts/vector.tcl

mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
set out [open $outFile w]
set cube [list {-1 0 0} {1 0 0} {0 -1 0} {0 1 0} {0 0 -1} {0 0 1}]

foreach quiet {0} {
    set posList [$sel get {x y z}]
    set radList [$sel get radius]
}

foreach pos $posList rad $radList {
    foreach d $cube {
	puts $out [vecAdd $pos [vecScale $rad $d]]
    }
    
    for {set i 0} {$i < $points} {incr i} {
	puts $out [vecAdd $pos [vecScale $rad [vecRandom]]]
    }
}

close $out
exit
