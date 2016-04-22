# This script will remove water from psf and pdf files.
# Use with: vmd -dispdev text -e removeWater.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "resname SIN"
#Input:
set psf pore1.2_all.psf
set pdb pore1.2_all.pdb
#Output:
set finalPsf pore1.2_all_types.psf

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb

set sel [atomselect top "($selText) and name \"SI.*\""]
set newType {}
foreach nb [$sel get numbonds] {
    lappend newType SI_${nb}
}
$sel set type $newType
$sel set charge 0.771000
$sel delete

set sel [atomselect top "($selText) and name \"N.*\""]
set newType {}
foreach nb [$sel get numbonds] {
    lappend newType N_${nb}
}
$sel set type $newType
$sel set charge -0.578359
$sel delete

set all [atomselect top all]
$all writepsf $finalPsf
puts "Types: [lsort -unique [$all get type]]"
$all delete
exit



