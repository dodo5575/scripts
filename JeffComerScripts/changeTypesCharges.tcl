# This script will change the types of atoms.
# Use with: vmd -dispdev text -e changeTypes.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set changeList {{"type O" OSI}}
set posCharge 2.4
set posType SI
set negCharge -1.2
set negType OSI
#Input:
set psf SiO2_r15z15.psf
set pdb SiO2_r15z15.pdb
#Output:
set finalPsf SiO2_square.psf

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

# Set the charges.
set sel [atomselect top "type $posType"]
$sel set charge $posCharge
$sel delete
set sel [atomselect top "type $negType"]
$sel set charge $negCharge
$sel delete

$all writepsf $finalPsf
puts "Final types: [lsort -unique [$all get type]]"
$all delete
exit



