# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "type \"N_.*\""
set selTextLeaving "type SI_4"
#Input:
set psf pore22.psf
set pdb pore22.pdb
#Output:
set finalPsf pore22_neutral.psf
set finalPdb pore22_neutral.pdb

# Load the molecule.
mol load psf $psf pdb $pdb
set all [atomselect top all]
set sel [atomselect top $selText]

# Get the initial charge.
set charge [measure sumweights $all weight charge]
puts "Initial charge: $charge"

# Distribute it among selText.
set qi [lindex [$sel get charge] 0]
set num [$sel num]
set q [expr -$charge/$num + $qi]
set err [expr ($q-$qi)/$qi]
puts "Shifting the charge of `$selText' from $qi to $q"
puts "Relative error: $err"
$sel set charge $q

# Write the intermediate result.
$all writepsf $finalPsf

# Clean up.
$sel delete
$all delete
mol delete top

# Reload to fix round off error.
mol load psf $finalPsf pdb $pdb
set all [atomselect top all]
set newCharge [measure sumweights $all weight charge]
puts "New charge: $newCharge"
set sel [atomselect top $selTextLeaving]
set q [expr [lindex [$sel get charge] 0] - $newCharge]
set index [lindex [$sel get index] 0]
$sel delete

# Shift.
puts "Shifting charge of atom $index by [expr -$newCharge]."
set sel [atomselect top "index $index"]
$sel set charge $q
set newCharge [measure sumweights $all weight charge]
puts "Final charge: $newCharge"

# Write the results.
$all writepdb $finalPdb
$all writepsf $finalPsf

$sel delete
$all delete
mol delete top
exit



