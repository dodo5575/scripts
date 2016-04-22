# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "segname ADNA BDNA and name P"
# Input:
set psf extra_water7226.psf
set coord extra_water7226.pdb

mol load psf $psf
mol addfile $coord
set sel [atomselect top $selText]
set cen [measure center $sel weight mass]
puts ""
puts "center: "
puts $cen

$sel delete
mol delete top
exit



