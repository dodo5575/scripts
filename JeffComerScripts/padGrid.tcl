# Input:
set inFile third_grisha_pad.dx
# Output:
set outFile third_grisha_in.dx

source vector.tcl
source gridForce.tcl

readDx grid0 $inFile
puts "Read $inFile, containing [llength $grid0(data)] potential values."

puts "Padding..."
padGrid grid0 grid4 {-2 -2 0}
puts "Initial size: $grid0(nx) $grid0(ny) $grid0(nz)"
puts "Final size: $grid4(nx) $grid4(ny) $grid4(nz)"

puts "Writing $outFile with [llength $grid4(data)] potential values."
writeDx grid4 $outFile
puts "Done."



