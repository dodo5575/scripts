# Cut a parabolic pore.
# use with: vmd -dispdev text -e cutParabola.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

#Input:
set pdb hcp40.pdb
set psf hcp40.psf
#Output:
set finalpdb pore.pdb
set finalpsf pore.psf
# Parameters:
set rmin 50

mol load pdb $pdb
set all [atomselect top all]

# Center the system.
set cen [measure center $all weight mass]
$all moveby [vecinvert $cen]
$all writepdb $finalpdb

# Check the system geometry.
set minmax [measure minmax $all]
set dx [expr [lindex $minmax 1 0] - [lindex $minmax 0 0]]
set dy [expr [lindex $minmax 1 1] - [lindex $minmax 0 1]]
set dz [expr [lindex $minmax 1 2] - [lindex $minmax 0 2]]
set zc [expr [lindex $minmax 0 2] + 0.5*$dz]
if {$dx < $dy} {
	set rmax 0.5*$dx*0.9
} else {
	set rmax 0.5*$dy*0.9
}

# Select the parabola.
set a [expr 4.*($rmax-$rmin)/($dz*$dz)]
set selText "x^2 + y^2 +z^2 < ($a*(z-$zc)^2 + $rmin)^2"
set bad [atomselect top $selText]
set atomList [lsort -unique [$bad get {segid resid name}]]

# Use psfgen to delete the atoms within the parabola.
package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $finalpdb

# Delete the selection.
foreach atom $atomList {
	delatom [lindex $atom 0] [lindex $atom 1] [lindex $atom 2]
}

writepsf $finalpsf
writepdb $finalpdb
exit





