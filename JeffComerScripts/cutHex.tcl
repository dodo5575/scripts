# This script will remove water from psf and pdf in a
# hexagonal prism along the z-axis with a vertex along
# the x-axis.
# It also suggests some basis vectors.
# use with: vmd -dispdev text -e myCutWaterHex.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>
package require psfgen 1.3

# The water is cut from the top and bottom at a distance of "z0" from the center.
set z0 25.0
set a 60.0
# This is the stuff that is removed.
set waterText "all"
set swapXY 1

# Input:
set psf ASiO2.psf
set pdb ASiO2_cen.pdb
# Output:
set finalpsf hex_r6h5.psf
set finalpdb hex_r6h5.pdb

# Load the molecule.
mol load psf $psf pdb $pdb

if {!$swapXY} {
    set x x
    set y y
} else {
    set x y
    set y x
}

# Find water outside of the hexagon.
set sqrt3 [expr sqrt(3.0)]
set r [expr $a/$sqrt3]
# Check the middle rectangle.
set check "($waterText) and ((abs(${x}) < 0.5*$r and abs(${y}) > 0.5*$sqrt3*$r) or"
# Check the lines forming the nonhorizontal sides.
set check "$check (${y} > $sqrt3*(${x}+$r) or ${y} < $sqrt3*(${x}-$r) or"
set check "$check ${y} > $sqrt3*($r-${x}) or ${y} < $sqrt3*(-${x}-$r)) or"
# Cut the top and bottom.
set check "$check abs(z) > $z0)"
set w [atomselect top $check]
foreach null {0} {set violators [lsort -unique [$w get {segname resid name}]]}

# Remove the offending water molecules.
puts "Deleting the offending water molecules..."
resetpsf
readpsf $psf
coordpdb $pdb
foreach waterMol $violators {
	delatom [lindex $waterMol 0] [lindex $waterMol 1] [lindex $waterMol 2]
}

writepsf $finalpsf
writepdb $finalpdb
exit



