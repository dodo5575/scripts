# Move the center of mass to the origin.
# use with: vmd -dispdev text -e center.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>
set selText "resname PCGL OLEO"
set wrapText "water or ions"
# Input:
set psf bilayer_dopc_1M.psf
set coor run2_dopc_wat.restart.coor
set xsc run2_dopc_wat.restart.xsc
# Output:
set finalpdb run2_dopc_wat.pdb

package require pbctools

mol load psf $psf
mol addfile $coor waitfor all

# Center the selection.
set all [atomselect top "all"]
set sel [atomselect top $selText]
set cen [measure center $sel weight mass]
$all moveby [vecinvert $cen]

# Wrap based on the new positions.
pbc readxst $xsc
pbc wrap -compound res -sel $wrapText -all -shiftcenterrel {-0.5 -0.5 -0.5}

$all writepdb $finalpdb
$all delete
$sel delete
exit
