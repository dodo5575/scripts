# Move the center of mass to the origin.
# use with: vmd -dispdev text -e center.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>
set selText "segname ADNA BDNA"

# Input:
set psf DNA_start.psf
set coor equilibrate-30.restart.coor
# Output:
set finalpdb double_30.pdb

mol load psf $psf
mol addfile $coor waitfor all
set all [atomselect top "all"]
set sel [atomselect top $selText]
set cen [measure center $sel weight mass]
$all moveby [vecinvert $cen]
$all writepdb $finalpdb
$all delete
$sel delete
exit



