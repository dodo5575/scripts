set pore 2.0
# Input:
set inFile0 pore_grisha_third_inside.dx
set inFile1 pore_grisha_third_surf.dx
set scale1 5.0
# Output:
set outFile third_grisha_init.dx

source vector.tcl
source gridForce.tcl

readDx grid0 $inFile0
puts "Read $inFile0, containing [llength $grid0(data)] potential values."
readDx grid1 $inFile1
puts "Read $inFile1, containing [llength $grid1(data)] potential values."

scalePotential grid1 $scale1
puts "Scaled potential from $inFile1 by $scale1."
puts "Adding grids..."
addGridsSame grid0 grid1 grid2

puts "Trimming grids..."
trimGrid grid2 grid3
puts "Initial size: $grid2(nx) $grid2(ny) $grid2(nz)"
puts "Final size: $grid3(nx) $grid3(ny) $grid3(nz)"

puts "Padding..."
padGrid grid3 grid4 {0 0 2}
puts "Initial size: $grid3(nx) $grid3(ny) $grid3(nz)"
puts "Final size: $grid4(nx) $grid4(ny) $grid4(nz)"

puts "Writing $outFile with [llength $grid4(data)] potential values."
writeDx grid4 $outFile
puts "Done."



