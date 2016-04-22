# Scale the system to the desired size.
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "all"
set scaleX [expr 1.0]
set scaleY [expr 1.0]
set scaleZ [expr 71.1161/71.5097]
#Input:
set psf pore_no_basepair_1M.psf
set coor output/pore_no_1M_eq1.restart.coor
#Output:
set finalPdb pore_no_1M_eq1.pdb

mol load psf $psf
mol addfile $coor
# Scale the positions.
set sel [atomselect top $selText]
foreach zero {0} {set pos [$sel get {x y z}]}
set newPos {}
foreach r $pos {
    foreach {x y z} $r {break}
    set x [expr {$scaleX*$x}]
    set y [expr {$scaleY*$y}]
    set z [expr {$scaleZ*$z}]
    lappend newPos [list $x $y $z]
}
$sel set {x y z} $newPos
$sel delete

# Write the results.
set all [atomselect top all]
$all writepdb $finalPdb
$all delete

mol delete top
exit



