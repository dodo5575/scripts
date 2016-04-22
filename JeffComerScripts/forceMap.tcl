set filePrefix pore1.6_open_4V_cyl
set dl 6.0
set x0 0.0

# Input:
set inFile ${filePrefix}.dx
# Output:
set outFile open_p1.6_cyl_force.txt

source vector.tcl
source gridForce.tcl

readDx grid $inFile
puts "Read $inFile, containing [llength $grid(data)] potential values."
puts "Size: $grid(nx) $grid(ny) $grid(nz)"

set origin $grid(origin)
set d [getGridSize grid]
set nx [expr int(ceil([lindex $d 0]/$dl))]
set ny [expr int(ceil([lindex $d 1]/$dl))]
set nz [expr int(ceil([lindex $d 2]/$dl))]

puts "Writing force map of [expr $ny*$nz] points to $outFile..."
set out [open $outFile w]
set r0 [list $x0 [lindex $origin 1] [lindex $origin 2]]
for {set iy 0} {$iy < $ny} {incr iy} {
    for {set iz 0} {$iz < $nz} {incr iz} {
	set r [vecAdd $r0 [list 0.0 $iy*$dl $iz*$dl]]
	set fx [interpolateForce grid $r 0]
	set fy [interpolateForce grid $r 1]
	set fz [interpolateForce grid $r 2]
	foreach {x y z} $r {break}

	puts $out "$x $y $z $fx $fy $fz"
    }
}
close $out
puts "Done."



