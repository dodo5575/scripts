set filePrefix phantom5-5_pore2.2
set dl 8.0
set y0 0.0

# Input:
set inFile ${filePrefix}.dx
# Output:
set outFile ${filePrefix}_force.txt

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
set r0 [list [lindex $origin 0] $y0 [lindex $origin 2]]
for {set ix 0} {$ix < $nx} {incr ix} {
    for {set iz 0} {$iz < $nz} {incr iz} {
	set r [vecAdd $r0 [list $ix*$dl 0.0 $iz*$dl]]
	set fx [interpolateForce grid $r 0]
	set fy [interpolateForce grid $r 1]
	set fz [interpolateForce grid $r 2]
	foreach {x y z} $r {break}

	puts $out "$x $y $z $fx $fy $fz"
    }
}
close $out
puts "Done."



