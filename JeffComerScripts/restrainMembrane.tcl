# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
# Spring constant in kcal/(mol A^2)
set beta 20
set selText "resname SIO2"
set surfText "not all"
# Input:
set psf trap2.2_KCl.psf
set pdb trap2.2_KCl.pdb
set coordPdb silicon_rest20.pdb
# Output:
set restFilePrefix restrain_silica

# Get the desired coordinates.
mol load pdb $coordPdb
set sel [atomselect top $selText]
foreach zero {0} {set pos [$sel get {x y z}]}
$sel delete
mol delete top

mol load psf $psf pdb $pdb
set selAll [atomselect top all]

# Set the spring constants to zero for all atoms.
$selAll set occupancy 0.0
$selAll set beta 0.0

# Select the membrane.
set selMem [atomselect top $selText]
# Set the positions.
$selMem set {x y z} $pos

# Select the surface.
set selSurf [atomselect top "(${selText}) and (${surfText})"]

# Set the spring constant for the membrane to this beta value.
$selMem set beta $beta
# Constrain the surface 10 times more than the bulk.
$selSurf set beta [expr 10.0*$beta]
# Write the constraint file.
$selAll writepdb ${restFilePrefix}${beta}.pdb

$selMem delete
$selSurf delete
$selAll delete
mol delete top
exit



