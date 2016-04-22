# externalForce.tcl
# Add an external electric field to a potential grid.
# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set lx 30
set ly 30
set lz 240
set diam0 16
set diam1 20
set dl 0.5

source vector.tcl
source gridForce.tcl

set s0 [expr $diam0*0.5]
set s1 [expr $diam1*0.5]
set sm [expr 0.5*($s0 + $s1)]
set ds [expr $s1 - $s0]
set ds3 [expr $ds*$ds*$ds]

proc getPotential {r} {
    global s0 s1 ds sm ds3
    foreach {x y z} $r {break}

    set s [expr sqrt($x*$x + $y*$y)]

    if {$s >= $s1} {return 1}
    if {$s <= $s0} {return 0}

    set sr [expr $s-$sm]
    return [expr -2.0/$ds3*$sr*$sr*$sr + 1.5/$ds*$sr + 0.5]
}

newGridDim g $dl $dl [expr ceil($lz/5)] $lx $ly $lz
set pot {}
for {set ix 0} {$ix < $g(nx)} {incr ix} {
    for {set iy 0} {$iy < $g(ny)} {incr iy} {
	for {set iz 0} {$iz < $g(nz)} {incr iz} {
	    set r [getPosition g [list $ix $iy $iz]]
	    lappend pot [getPotential $r]
	}
    }
}
set g(data) $pot

puts "origin: $g(origin)"
puts "number of nodes: [llength $g(data)]"

writeDx g cylinder${diam0}-${diam1}_sys$lx-$ly.dx
exit

