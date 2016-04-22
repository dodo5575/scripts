# potentialProfileAverage.tcl
# Average the potential profile in slices in the xy-plane, obtaining V(z).
# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@illinois.edu>

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

if {$argc != 2} {
    puts "This script requires two arguments: "
    puts "  radius over which to average"
    puts "  file name prefix"
    puts ""
    exit
}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}

set radius [lindex $argv 0]
set filePrefix [trimExtension [lindex $argv 1]]
# Input:
set inFile ${filePrefix}.dx
# Output:
set outFile ${filePrefix}_cylz.dat

set radius2 [expr $radius*$radius]
readDx grid $inFile
puts "Read $inFile, containing [llength $grid(data)] potential values."
puts "Size: $grid(nx) $grid(ny) $grid(nz)"

puts "Writing potential profile of $grid(nz) points to $outFile..."
set out [open $outFile w]
for {set iz 0} {$iz < $grid(nz)} {incr iz} {
    set r [getPosition grid [list 0 0 $iz]]
    set z [lindex $r 2]
    set pot 0.0

    # Compute the average potential in this slice.
    set n 0
    for {set iy 0} {$iy < $grid(ny)} {incr iy} {
	for {set ix 0} {$ix < $grid(nx)} {incr ix} {
	    set r [getPosition grid [list $ix $iy $iz]]
	    foreach {x y z} $r {break}
	    if {$x*$x + $y*$y < $radius2} {
		set i [expr $iz + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix]
		set pot [expr $pot + [lindex $grid(data) $i]]
		incr n
	    }
	}
    }
   
    set pot [expr $pot/$n]
    puts $out "$z $pot"
}
close $out
puts "Done."



