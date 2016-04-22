# Apply a constant force to the water.
# to use: vmd -dispdev text -e pushWater.tcl

# Push the water with this force in kcal/(mol*A).
# You cannot record a force smaller than 0.00001 in the PDB.
set force 0.0005
set selText "water or ions"
# Input:
set psf system-ions.psf
set pdb system-ions.pdb
# Output:
set forceFile waterForce_${force}.pdb

mol load psf $psf pdb $pdb
set selAll [atomselect top all]
$selAll set occupancy 0.0

set sel [atomselect top $selText]
$sel set x 0.
$sel set y 0.
$sel set z [expr -$force*100.]
# Use the smallest scaling factor possible.
$sel set occupancy 0.01

$selAll writepdb $forceFile
exit



