# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e markDna.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
# Spring constant in kcal/(mol A^2)
set selText "segname HAIR"
# Input:
set psf phantom_p2.0_bsc.psf
set pdb phantom_p2.0_bsc.pdb
# Output:
set restFilePrefix mark_dna_charge

mol load psf $psf pdb $pdb
set selAll [atomselect top all]

# Set beta to zero for all atoms.
$selAll set beta 0.0
$selAll set occupancy 0.0

# Select the atoms to which langevin is applied.
set sel [atomselect top $selText]
foreach zero {0} {set q [$sel get charge]}

# Set the beta column.
$sel set beta 1.0
$sel set occupancy $q
# Write the file.
$selAll writepdb ${restFilePrefix}.pdb
puts "Marked [$sel num] atoms with 1.0 in B and the charge in O."

$sel delete
$selAll delete
mol delete top
exit



