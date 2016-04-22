# This script will remove water from psf and pdf in a
# hexagonal prism along the z-axis with a vertex along
# the x-axis.
# It also suggests some basis vectors.
# use with: vmd -dispdev text -e myCutWaterHex.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>
package require psfgen 1.3

# The water is cut from the top and bottom at a distance of "z0" from the center.
set z0 [expr 50+1000]
# The radius of the water hexagon is reduced by "radiusMargin"
# from the pore hexagon. The distance is in angstroms.
set radiusMargin 0.0
# This is the stuff that is removed.
set waterText "water or ions"
# This selection forms the basis for the hexagon.
set selText "resname SIN"
# The angle at which we find a flat face of the hexagon.
set angle 90

# Input:
set psf cont_dna_sol.psf
set pdb cont_dna_sol.pdb
# Output:
set finalpsf cont_dna_hex.psf
set finalpdb cont_dna_hex.pdb
set basistext cellbasis.txt

# Load the molecule.
mol load psf $psf pdb $pdb

# Put the origin at the center of the pore.
set sel [atomselect top $selText]
set cen [measure center $sel weight mass]
set all [atomselect top "all"]
$all moveby [vecinvert $cen]
$all move [transaxis z $angle]

# Find the system dimensions.
set minmax [measure minmax $sel]
set size [vecsub [lindex $minmax 1] [lindex $minmax 0]]
foreach {size_x size_y size_z} $size {break}
# This is the hexagon's radius.
set rad [expr 0.5*$size_y]
set r [expr $rad - $radiusMargin]

# Find water outside of the hexagon.
set sqrt3 [expr sqrt(3.0)]
# Check the middle rectangle.
set check "($waterText) and ((abs(y) < 0.5*$r and abs(x) > 0.5*$sqrt3*$r) or"
# Check the lines forming the nonhorizontal sides.
set check "$check (x > $sqrt3*(y+$r) or x < $sqrt3*(y-$r) or"
set check "$check x > $sqrt3*($r-y) or x < $sqrt3*(-y-$r)) or"
# Cut the top and bottom.
set check "$check abs(z) > $z0)"
set w [atomselect top $check]
foreach null {0} {set violators [lsort -unique [$w get {segname resid}]]}

# Announce the suggested cell basis vectors.
puts ""
puts "Suggested cell basis vectors: (one unit cell too small!)"
set cellBasisVector1 [list [expr 1.5*$rad] [expr 0.5*$sqrt3*$rad] 0.0]
set cellBasisVector2 [list 0.0 [expr $sqrt3*$rad] 0.0]
set cellBasisVector3 [list 0.0 0.0 $size_z]
puts ""

# Open the file and write the suggested cell basis vectors.
set out [open $basistext w]
puts $out "These vectors are one unit cell too small."
puts $out [format "cellBasisVector1 %s" $cellBasisVector1]
puts $out [format "cellBasisVector2 %s" $cellBasisVector2]
puts $out [format "cellBasisVector3 %s" $cellBasisVector3]
close $out

# Remove the offending water molecules.
puts "Deleting the offending water molecules..."
resetpsf
readpsf $psf
coordpdb $pdb
foreach waterMol $violators {
	delatom [lindex $waterMol 0] [lindex $waterMol 1]
}

writepsf $finalpsf
writepdb $finalpdb
exit



