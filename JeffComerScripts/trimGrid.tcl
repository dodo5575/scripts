set filePrefix trap2.0_init_third
# Input:
set inFile $filePrefix.dx
set outFile ${filePrefix}_trim.dx

source vector.tcl
source gridForce.tcl

readDx grid $inFile
puts "Read $inFile, containing [llength $grid(data)] potential values."
puts "Size: $grid(nx) $grid(ny) $grid(nz)"

trimGrid grid new
puts "Writing $outFile with [llength $new(data)] potential values."
puts "Size: $new(nx) $new(ny) $new(nz)"
writeDx new $outFile
puts "Done."



