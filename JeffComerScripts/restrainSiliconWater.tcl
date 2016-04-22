# Add harmonic constraints to silicon nitride.
# jcomer2@uiuc.edu
#set sys [lindex $argv 0]
#set water [lindex $argv 1]

set water 110413
set sys s-dna
#set water 91924
#set sys b-dna
# Parameters:
# Spring constant in kcal/(mol A^2)
set betaList {1}
set selText "resname SIN"
set surfText "(name \"SI.*\" and numbonds<=3) \
or (name \"N.*\" and numbonds<=2)"
# Input:
set psf pore6_${sys}_water${water}.psf
set pdb pore6_${sys}_water${water}.pdb
set coordPsf pore6_${sys}_conc.psf
set coordPdb pore6_${sys}_conc.pdb
# Output:
set restFilePrefix restrain_${sys}_sin_water${water}

mol load psf $coordPsf pdb $coordPdb
set sel [atomselect top $selText]
foreach quiet {0} { set pos [$sel get {x y z}] }
$sel delete
mol delete top

mol load psf $psf pdb $pdb
set selAll [atomselect top all]

# Set the spring constants to zero for all atoms.
$selAll set occupancy 0.0
$selAll set beta 0.0

# Select the silicon nitride.
set selSiN [atomselect top $selText]
# Set the coordinates.
$selSiN set {x y z} $pos

# Select the surface.
set selSurf [atomselect top "(${selText}) and (${surfText})"]

foreach beta $betaList {
	# Set the spring constant for SiN to this beta value.
	$selSiN set beta $beta
	# Constrain the surface 10 times more than the bulk.
	$selSurf set beta [expr 10.0*$beta]
	# Write the constraint file.
	$selAll writepdb ${restFilePrefix}.pdb
}
$selSiN delete
$selSurf delete
$selAll delete
mol delete top
exit
