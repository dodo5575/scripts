# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e langevinWater.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
# Spring constant in kcal/(mol A^2)
set betaList {0.2}
set selText "water and name OH2"
# Input:
set psf phtm_tail_all.psf
set pdb phtm_tail_all.pdb
# Output:
set restFilePrefix langevin

mol load psf $psf pdb $pdb
set selAll [atomselect top all]

# Set the spring constants to zero for all atoms.
$selAll set occupancy 0.0
$selAll set beta 0.0

# Select the atoms to which langevin is applied.
set sel [atomselect top $selText]

foreach beta $betaList {
	# Set the Langevin coupling.
	$sel set beta $beta
	# Write the constraint file.
	$selAll writepdb ${restFilePrefix}_${beta}.pdb
}
$selSel delete
$selAll delete
mol delete top
exit



