# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set nRemove 500
# Parameters:
set selText "type SI and (x^2+y^2<(10 + 0.36397*abs(z))^2 or abs(z) > 55)"
#Input:
set psf SiO2_square.psf
set pdb SiO2_square.pdb
#Output:
set finalPsf silica_square.psf
set finalPdb silica_square.pdb

# Load the system.
mol load psf $psf pdb $pdb
set all [atomselect top all]

# Get the charge.
#set q [measure sumweights $all weight charge]
#set nRemove [expr int(floor($q+0.5))]

set sel [atomselect top $selText]
set nSelect [$sel num]
puts ""
puts "Will delete $nRemove atoms from a set of $nSelect ions with ($selText)"

# Determine the total charge.
set charge [measure sumweights $all weight charge]
puts "The total charge is $charge"
$all delete

foreach zero {0} {
    set posList [$sel get {x y z}]
    set atomList [$sel get {segname resid name}]
}
$sel delete

# Get the sorting weight.
set weightList {}
foreach pos $posList atom $atomList {
    foreach {x y z} $pos {break}
    set weight [expr $x*$x + $y*$y]
    lappend weightList [concat $weight $atom]
}
# Sort the atoms.
foreach zero {0} {
    set sortList [lsort -real -index 0 $weightList]
    set remove [lrange $sortList 0 [expr $nRemove-1]]
}

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $pdb

# Delete the selection.
foreach atom $remove {
    delatom [lindex $atom 1] [lindex $atom 2] [lindex $atom 3]
}

writepsf $finalPsf
writepdb $finalPdb
puts ""
puts "[llength $remove] atoms were deleted."

# Determine the total charge of the final system.
mol delete top
mol load psf $finalPsf pdb $finalPdb
set all [atomselect top all]
set chargeNew [measure sumweights $all weight charge]
puts "Initial charge $charge"
puts "Final charge $chargeNew"
$all delete
mol delete top
exit
