# potentialScale.tcl
# Scale by a factor.
# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 1} {
	puts "This script requires two arguments: "
	puts "  file name"
        puts ""
	exit
}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}

# Parameters:
set prefix [trimExtension [lindex $argv 0]]
set outWorld 6

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

set inFile ${prefix}.dx

puts "Reading $inFile"
readDx grid $inFile

puts "Running"
set slope {}
set inter {}
for {set ix 0} {$ix < $grid(nx)} {incr ix} {
    for {set iy 0} {$iy < $grid(ny)} {incr iy} {
	set iz0 [expr {$grid(nz)-1-$outWorld}]
	set iz1 $outWorld
	set j0 [expr {$iz0 + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix}]
	set j1 [expr {$iz1 + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix}]
	set pot0 [lindex $grid(data) $j0]
	set pot1 [lindex $grid(data) $j1]

	set slope [expr {($pot1 - $pot0)/(2.0*$outWorld)}]
	
	for {set i 1} {$i <= $outWorld} {incr i} {
	    set iz [expr {$i + $grid(nz)-$outWorld-1}]
	    set j [expr {$iz + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix}]
	    set pot [expr {$pot0 + $slope*$i}]

	    lset grid(data) $j $pot
	}

	for {set i 0} {$i < $outWorld-1} {incr i} {
	    set iz $i
	    set j [expr {$iz + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix}]
	    set pot [expr {$pot0 + $slope*($i+$outWorld)}]

	    lset grid(data) $j $pot
	}
	
    }
}

puts "Writing"
writeDx grid ${prefix}_edges.dx



