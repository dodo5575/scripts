set all [atomselect top "water"]
set minmax [measure minmax $all]
set dx [expr [lindex $minmax 1 0]-[lindex $minmax 0 0]]
set dy [expr [lindex $minmax 1 1]-[lindex $minmax 0 1]]
set dz [expr [lindex $minmax 1 2]-[lindex $minmax 0 2]]
molinfo top set alpha 90
molinfo top set beta 90
molinfo top set gamma 90
molinfo top set a $dx
molinfo top set b $dy
molinfo top set c $dz



