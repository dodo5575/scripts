# potentialProfileAverage.tcl
# Average the potential profile in slices in the xy-plane, obtaining V(z).
# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 2} {
    puts "This script requires two arguments: "
    puts "  file name"
    puts "  factor"
    puts ""
    exit
}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}

set filePrefix [trimExtension [lindex $argv 0]]
# Input:
set inFile [lindex $argv 0]
set factor [lindex $argv 1]
# Output:
set outFile ${filePrefix}_z.$factor.dat

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

readDx grid $inFile
puts "Read $inFile, containing [llength $grid(data)] potential values."
puts "Size: $grid(nx) $grid(ny) $grid(nz)"

set ix [expr int(floor($factor*$grid(nx) + 0.5))]
set iy [expr int(floor($factor*$grid(ny) + 0.5))]

puts "Writing potential profile of $grid(nz) points to $outFile..."
set out [open $outFile w]
for {set iz 0} {$iz < $grid(nz)} {incr iz} {
    set r [getPosition grid [list 0 0 $iz]]
    set z [lindex $r 2]
    
    set i [expr $iz + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix]
    set pot [lindex $grid(data) $i]

    puts $out "$z $pot"
}
close $out
puts "Done."
