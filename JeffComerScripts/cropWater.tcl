# This script will remove water from psf and pdf files.
# Use with: vmd -dispdev text -e cropWater.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "water or ions"
set margin 25

#Input:
set psf  1p2y.psf
set pdb  1p2y.pdb
set bincoords eq4.restart.coor
#Output:
set finalpsf 1p2y_reduced.psf
set finalpdb 1p2y_reduced.pdb
set temppdb tmp.pdb

mol load psf $psf pdb $pdb
mol addfile $bincoords waitfor all

# Write a pdb with coordinates from "bincoords".
set all [atomselect top all]
$all writepdb $temppdb

# Obtain the new water box size.
set other [atomselect top "not ($selText)"]
set minmax [measure minmax $other]
set x0 [expr [lindex $minmax 0 0] - $margin]
set y0 [expr [lindex $minmax 0 1] - $margin]
set z0 [expr [lindex $minmax 0 2] - $margin]
set x1 [expr [lindex $minmax 1 0] + $margin]
set y1 [expr [lindex $minmax 1 1] + $margin]
set z1 [expr [lindex $minmax 1 2] + $margin]

set violatorText "($selText) and (x < $x0 or x > $x1 or y < $y0 or y > $y1 \
or z < $z0 or z > $z1)" 
set sel [atomselect top $violatorText]
set atomList [lsort -unique [$sel get {segname resid}]]

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $temppdb

# Delete the selection.
foreach atom $atomList {
	delatom [lindex $atom 0] [lindex $atom 1]
}

writepsf $finalpsf
writepdb $finalpdb

puts "Dimensions of new water box: "
puts "[expr $x1-$x0] [expr $y1-$y0] [expr $z1-$z0]"
exit



