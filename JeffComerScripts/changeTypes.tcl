# This script will change the types of atoms.
# Use with: vmd -dispdev text -e changeTypes.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set changeList {{"type OT" OW} {"type HT" HW} {"name C5'" CI}}
#Input:
set psf pore2.0_all.psf
set pdb pore2.0_all.pdb
#Output:
set finalPsf pore2.0_CI.psf

mol load psf $psf pdb $pdb
set all [atomselect top all]
puts "Initial types: [lsort -unique [$all get type]]"

# Change the types.
foreach pair $changeList {
    set inText [lindex $pair 0]
    set outType [lindex $pair 1]
    set sel [atomselect top $inText]
    $sel set type $outType
    $sel delete
}

$all writepsf $finalPsf
puts "Final types: [lsort -unique [$all get type]]"
$all delete
exit



