# Author: Jeff Comer <jcomer2@illinois.edu>

# Input:
set inGrid diffusion70_pot.dx
set inFile grid_energy_low.log.dat
# Output:
set outFile phantom_low.dx

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

readDx grid $inGrid
puts "Read $inGrid, containing [llength $grid(data)] potential values."
puts "Size: $grid(nx) $grid(ny) $grid(nz)"

puts "Loading the data from $inFile."
set data {}
set in [open $inFile r]
while {[gets $in line] >= 0} {
    if {[string length $line] < 1} { continue }
    lappend data $line
}
close $in
puts "Done."

set grid(data) $data

puts "Wrote the result $outFile."
writeDx grid $outFile

