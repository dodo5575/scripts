# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "all"
# Input:
set psf cg_ions_0.4M_big.psf
set coord output/big_bulk_0.4M_eq0.restart.coor

mol load psf $psf
mol addfile $coord
set sel [atomselect top $selText]
set minmax [measure minmax $sel]
set dx [expr [lindex $minmax 1 0]-[lindex $minmax 0 0]]
set dy [expr [lindex $minmax 1 1]-[lindex $minmax 0 1]]
set dz [expr [lindex $minmax 1 2]-[lindex $minmax 0 2]]
puts ""
puts "System size: "
puts "$dx $dy $dz"

$sel delete
mol delete top
exit



