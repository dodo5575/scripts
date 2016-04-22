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
set outFile ${filePrefix}_y.$factor.dat

source vector.tcl
source gridForce.tcl

readDx grid $inFile
puts "Read $inFile, containing [llength $grid(data)] potential values."
puts "Size: $grid(nx) $grid(ny) $grid(nz)"

set iz [expr int(floor($factor*$grid(nz) + 0.5))]
set ix [expr int(floor($factor*$grid(nx) + 0.5))]

puts "Writing potential profile of $grid(nz) points to $outFile..."
set out [open $outFile w]
for {set iy 0} {$iy < $grid(ny)} {incr iy} {
    set r [getPosition grid [list 0 $iy 0]]
    set y [lindex $r 1]
    
    set i [expr $iz + $grid(nz)*$iy + $grid(nz)*$grid(ny)*$ix]
    set pot [lindex $grid(data) $i]

    puts $out "$y $pot"
}
close $out
puts "Done."
