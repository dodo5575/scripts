
puts "Usage: alignCenter {inpdbPrefix outpdbPrefix}" 

proc alignCenter {inpdbPrefix outpdbPrefix} {

set in [mol new $inpdbPrefix.pdb]

set all [atomselect top all]
set water [atomselect top water]

set box [measure minmax $water] 
set boxMinX [lindex $box 0 0]
set boxMinY [lindex $box 0 1]
set boxMinZ [lindex $box 0 2]
set boxMaxX [lindex $box 1 0]
set boxMaxY [lindex $box 1 1]
set boxMaxZ [lindex $box 1 2]

set inDNAselText "nucleic and (x > $boxMinX and y > $boxMinY and z > $boxMinZ) and (x < $boxMaxX and y < $boxMaxY and z < $boxMaxZ)"

set inDNA [atomselect top $inDNAselText]

set minMax [measure minmax $inDNA] 
set MinX [lindex $minMax 0 0]
set MinY [lindex $minMax 0 1]
set MinZ [lindex $minMax 0 2]
set MaxX [lindex $minMax 1 0]
set MaxY [lindex $minMax 1 1]
set MaxZ [lindex $minMax 1 2]


set X [expr ($MaxX + $MinX) / 2 * (-1)]
set Y [expr ($MaxY + $MinY) / 2 * (-1)]
set Z [expr ($MaxZ + $MinZ) / 2 * (-1)]


set vector [list $X $Y $Z]
puts $vector

$all moveby $vector

$all writepdb $outpdbPrefix.pdb

mol delete $in

}
