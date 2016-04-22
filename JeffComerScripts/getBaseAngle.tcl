# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set axis {0 0 1}
set selText "segname ADNA"
# Input:
set psf pore_at_basepair.psf
set coor pore_at_basepair.pdb

mol load psf $psf
mol addfile $coor
set sel [atomselect top $selText]

set resName [lindex [$sel get resname] 0]

set mole top
# Get the hexagon ring basis.
if {[string equal ADE $resName] || [string equal GUA $resName]} {
    set selX0 [atomselect $mole "($selText) and name C4"]
    set selX1 [atomselect $mole "($selText) and name N1"]
    set selY0 [atomselect $mole "($selText) and name C4"]
    set selY1 [atomselect $mole "($selText) and name C6"]
} else {
    set selX0 [atomselect $mole "($selText) and name C6"]
    set selX1 [atomselect $mole "($selText) and name N3"]
    set selY0 [atomselect $mole "($selText) and name C6"]
    set selY1 [atomselect $mole "($selText) and name C4"]
}

set rX0 [lindex [$selX0 get {x y z}] 0]
set rX1 [lindex [$selX1 get {x y z}] 0]
set rY0 [lindex [$selY0 get {x y z}] 0]
set rY1 [lindex [$selY1 get {x y z}] 0]

$selX0 delete
$selX1 delete
$selY0 delete
$selY1 delete

set ex [vecsub $rX1 $rX0]
set ex [vecscale [expr 1.0/[veclength $ex]] $ex]
set ey [vecsub $rY1 $rY0]
set ey [vecsub $ey [vecscale [vecdot $ey $ex] $ex]]
set ey [vecscale [expr 1.0/[veclength $ey]] $ey]
set ez [veccross $ex $ey]

set pi [expr {4.0*atan(1.0)}]
set comp [vecdot $ez $axis]
set angle [expr {180.0/$pi*acos($comp)}]
puts ""
puts "normal to base plane: $ez"
puts "angle between normal to the base plane and axis: $angle"


$sel delete
mol delete top
exit



