# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e markDna.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
# Spring constant in kcal/(mol A^2)
set selText1 "segname HAIR and resid 103 104"
set selText2 "resname SIN"
set selText3 "water or ions"

# Input:
set psf pore+dna_E4V.psf
set pdb steady_4V2.pdb
# Output:
set restFilePrefix mark_groups

mol load psf $psf pdb $pdb
set selAll [atomselect top all]

# Set beta to zero for all atoms.
$selAll set beta 0.0

# Set the occupancy to the charge.
$selAll set occupancy [$selAll get charge]

# Select the atoms to which langevin is applied.
set sel [atomselect top $selText1]
# Set the beta column.
$sel set beta 1.0
puts "Marked [$sel num] atoms with 1.0."
$sel delete

# Select the atoms to which langevin is applied.
set sel [atomselect top $selText2]
# Set the beta column.
$sel set beta 2.0
puts "Marked [$sel num] atoms with 2.0."
$sel delete

# Select the atoms to which langevin is applied.
set sel [atomselect top $selText3]
# Set the beta column.
$sel set beta 3.0
puts "Marked [$sel num] atoms with 3.0."
$sel delete

# Write the file.
$selAll writepdb ${restFilePrefix}.pdb
$selAll delete
mol delete top
exit



