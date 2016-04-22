# Fix inner atoms.
# to use: vmd -dispdev text -e fixInner.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
# Spring constant in kcal/(mol A^2)
set betaList {1.0}
set selText "abs(z) < 43 and x^2+y^2 > (17 + 2*(28-17)/102*abs(z))^2"
# Input:
set psf silica_shell.psf
set pdb silica_shell.pdb
# Output:
set restFile silica_shell_ready1.pdb

mol load psf $psf pdb $pdb
set selAll [atomselect top all]

# Set the spring constants to zero for all atoms.
$selAll set occupancy 0.0
$selAll set beta 0.0

# Select the surface.
set selSurf [atomselect top $selText]

foreach beta $betaList {
	# Set occupancy and beta.
	$selSurf set occupancy $beta
	$selSurf set beta $beta
	$selAll writepdb $restFile
}
$selSurf delete
$selAll delete
mol delete top
exit



