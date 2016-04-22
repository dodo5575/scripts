# Use with: vmd -dispdev text -e shuffleIons.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set waterText "name OH2 and (abs(z)>18 or x^2+y^2<10^2)"
set ionText "name POT CLA"
#Input:
set psf pore_at_basepair.psf
set pdb pore_at_eq1.pdb
#Output:
set outName pore_at_basepair1

proc swapResiduePos {seg0 res0 seg1 res1} {
    set s0 [atomselect top "segname $seg0 and resid $res0"]
    set s1 [atomselect top "segname $seg1 and resid $res1"]

    if {[$s0 num] <= 0 || [$s1 num] <= 0} {
	puts "Warning! No atoms in ${seg0}:${res0} or ${seg1}:${res1}."
	return 0
    }

    set r0 [measure center $s0 weight mass]
    set r1 [measure center $s1 weight mass]

    $s0 moveby [vecsub $r1 $r0]
    $s1 moveby [vecsub $r0 $r1]

    $s0 delete
    $s1 delete
    return 1
}

# Load the system.
mol load psf $psf pdb $pdb
set all [atomselect top all]

set waterSel [atomselect top $waterText]
set ionSel [atomselect top $ionText]

foreach zero {0} {
    set waterResList [$waterSel get {segname resid}]
    set ionResList [$ionSel get {segname resid}]
}
$waterSel delete
$ionSel delete
set n [llength $ionResList]

# Sort the water randomly.
set sortList {}
foreach res $waterResList {
    lappend sortList [concat $res [expr rand()]]
}
foreach zero {0} {
    set sortList [lsort -real -index 2 $sortList]
    set sortList [lrange $sortList 0 [expr $n-1]]
}

foreach ion $ionResList water $sortList {
    swapResiduePos [lindex $ion 0] [lindex $ion 1] [lindex $water 0] [lindex $water 1]
}

puts "Scrambled $n ions."
$all writepdb $outName.pdb

$all delete
mol delete top
exit


