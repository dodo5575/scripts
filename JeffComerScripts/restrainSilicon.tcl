# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# jcomer2@uiuc.edu

set size 12
# Parameters:
# Spring constant in kcal/(mol A^2)
set betaList {1}
set selText "resname SIN"
set surfText "(name \"SI.*\" and numbonds<=3) \
or (name \"N.*\" and numbonds<=2)"
# Input:
set psf dna_water67700_mem10.psf
set pdb dna_water67700_mem10.pdb
# Output:
set restFilePrefix restrain_sin

mol load psf $psf pdb $pdb
set selAll [atomselect top all]

# Set the spring constants to zero for all atoms.
$selAll set occupancy 0.0
$selAll set beta 0.0

# Select the silicon nitride.
set selSiN [atomselect top $selText]

# Select the surface.
set selSurf [atomselect top "(${selText}) and (${surfText})"]

foreach beta $betaList {
	# Set the spring constant for SiN to this beta value.
	$selSiN set beta $beta
	# Constrain the surface 10 times more than the bulk.
	$selSurf set beta [expr 10.0*$beta]
	# Write the constraint file.
	$selAll writepdb ${restFilePrefix}${beta}.pdb
}
$selSiN delete
$selSurf delete
$selAll delete
mol delete top
exit
