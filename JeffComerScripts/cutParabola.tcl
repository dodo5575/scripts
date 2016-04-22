# Cut a parabolic pore.
# use with: vmd -dispdev text -e cutParabola.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

#Input:
set pdb  bigHexagon-cg.pdb
#Output:
set finalpdb bigPore-cg.pdb
# Parameters:
set rmin 30

mol load pdb $pdb
set all [atomselect top all]

# Center the system.
set cen [measure center $all weight mass]
$all moveby [vecinvert $cen]

# Check the system geometry.
set minmax [measure minmax $all]
set dx [expr [lindex $minmax 1 0] - [lindex $minmax 0 0]]
set dz [expr [lindex $minmax 1 2] - [lindex $minmax 0 2]]
set zc [expr [lindex $minmax 0 2] + 0.5*$dz]
set rmax 0.5*$dx*0.8

set a [expr 4.*($rmax-$rmin)/($dz*$dz)]
set selText "x^2 + y^2 +z^2 > ($a*(z-$zc)^2 + $rmin)^2"
set good [atomselect top $selText]

$good writepdb $finalpdb
exit






