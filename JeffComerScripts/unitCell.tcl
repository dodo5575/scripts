# This script renders the Si3N4 unit cell.
# use with: vmd -e unitCell.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

mol load pdb unitCell_alpha-unique.pdb

graphics top delete all
set sel [atomselect top all]
set cen [measure center $sel]

set amag 7.595
set bmag 7.595
set cmag 2.902
set aunit {1. 0. 0.}
set bunit [list 0.5 [expr 0.5*sqrt(3.)] 0.]
set cunit {0. 0. 1.}

molinfo top set a $amag
molinfo top set b $bmag
molinfo top set c $cmag
molinfo top set alpha 90.
molinfo top set beta 90.
molinfo top set gamma 60.

set a [vecscale $aunit $amag]
set b [vecscale $bunit $bmag]
set c [vecscale $cunit $cmag]
set ap [vecadd $b $c]
set bp [vecadd $c $a]
set cp [vecadd $a $b]
set op [vecadd $a $b $c]

set rc [vecscale $op 0.5]
set vertList {}
foreach r [list [veczero] $a $b $c $ap $bp $cp $op] {
	lappend vertList $r
	#lappend vertList [vecsub $r $rc]
}

$sel moveby [vecsub $rc $cen]

set edgeList {0 1 1 5 3 5 0 3 4 7 2 4 6 7 2 6 3 4 5 7 1 6 0 2}
graphics top color white
foreach {i j} $edgeList {
	graphics top line [lindex $vertList $i] [lindex $vertList $j]
}




