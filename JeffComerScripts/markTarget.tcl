# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e markDna.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
# Spring constant in kcal/(mol A^2)
set beta 1.0
set selText "segname ADNA BDNA PRTA PRTB"
# Input:
set psf bam_phantom_pull.psf
set pdb bam_phantom_pull.pdb
# Output:
set restFilePrefix mark_target

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
