# Apply a random rotation matrix with a uniform distribution
# on a sphere to the selection.
# Use with: vmd -dispdev text -e randomRotation.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "segname HAIR"
# The number of random rotations to produce.
set n 5
#Input:
set psf conform3.psf
set pdb conform3.pdb
#Output:
set prefix conform3_rotate

set pi [expr 4.*atan(1.)]

# Obtain the {segid resid name} for the selection.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
set all [atomselect top all]

# Find the center of mass of the selection.
set cen [measure center $sel weight mass]

for {set i 0} {$i < $n} {incr i} {
	$sel moveby [vecinvert $cen]	
	
	# Create a random rotation matrix, uniform on a sphere.
	set psi [expr 2.0*rand()*$pi]
	set phi [expr acos(2.0*rand()-1.0)]
	set theta [expr 2.0*rand()*$pi]
	set basis [transaxis z $psi rad]
	set basis [transmult [transaxis x $phi rad] $basis]
	set basis [transmult [transaxis z $theta rad] $basis]
	
	# Rotate and move back.
	$sel move $basis
	$sel moveby $cen
	
	set check [atomselect top "not ($selText) and within 2.0 of ($selText)"]
	if {[$check num] > 0} {puts "Warning: Atoms too close!"}
	
	# Write the file.
	$all writepdb "${prefix}${i}.pdb"
}

exit



