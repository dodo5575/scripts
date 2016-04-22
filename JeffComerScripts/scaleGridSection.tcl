set inFile pore1.2_eq_phtm1.6-1.8.dx
set outFile pore1.2_surf_phtm1.6-1.8.dx
set scaleFactor 5.0
set surfZ 44.0
set switchDist 3.0

source vector.tcl
source gridForce.tcl

readDx grid $inFile
puts "Read $inFile, containing [llength $grid(data)] potential values."

puts "Modifying grid..."
copyGridDim new grid
set new(data) {}
set n 0
foreach pot $grid(data) {
    set r [getPosition new [latticePosition new $n]]
    set z [expr abs([lindex $r 2])]
    if {abs($z) < $surfZ} {
	# Do not scale.
	lappend new(data) $pot
    } elseif {abs($z) < $surfZ + $switchDist} {
	# Linearly interpolate the scale over switchDist.
	set s [expr ($scaleFactor-1.0)/$switchDist*($z-$surfZ) + 1.0]
	set p [expr $pot*$s]
	lappend new(data) $p
    } else {
	# Apply a constant scale.
	lappend new(data) [expr $pot*$scaleFactor]
    }
    incr n
}

puts "Writing $outFile with [llength $newGrid(data)] potential values."
writeDx newGrid $outFile
puts "Done."





