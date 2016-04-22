# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e markDna.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
# Spring constant in kcal/(mol A^2)
set beta 1.0
set selText "segname ADNA and not name H4' H2'' H2' H5'' H5' H1' H3' C4' O4' C1' C2' P O1P O2P O5' C5' C3' O3' H5T H3T"
# Input:
set psf polyC10AC14.psf
set pdb polyC10AC14.pdb
# Output:
set restFilePrefix polyC10AC14_fix.pdb

mol load psf $psf pdb $pdb
set selAll [atomselect top all]

# Set beta to zero for all atoms.
$selAll set beta 0.0
$selAll set occupancy 0.0

# Select the atoms to which langevin is applied.
set sel [atomselect top $selText]

# Set the beta column.
$sel set beta $beta
$sel set occupancy $beta
# Write the file.
$selAll writepdb ${restFilePrefix}.pdb
puts "Marked [$sel num] atoms."

$sel delete
$selAll delete
mol delete top
exit



