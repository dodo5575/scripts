# Input:
set inCoords pore+dna_coords_all.txt
set prefix phantom_pore2.2_high
set inFile ${prefix}.dx
# Output:
set outFile ${prefix}_force_all.txt

source vector.tcl
source gridForce.tcl

proc readPoints {fileName} {
    set in [open $fileName r]
    set ret {}

    foreach line [split [read $in] \n] {
	if {[string equal [string index $line 0] "\#"]} {continue}
	set tok [concat $line]
	if {[llength $tok] < 3} {continue}
	lappend ret [lrange $tok 0 2]
    }
    return $ret
}

readDx grid $inFile
puts "Read $inFile, containing [llength $grid(data)] potential values."
puts "Size: $grid(nx) $grid(ny) $grid(nz)"

set point [readPoints $inCoords]
puts "Read $inCoords, containing [llength $point] points."

puts "Writing the forces to $outFile..."
set out [open $outFile w]
foreach r $point {
    set fx [interpolateForce grid $r 0]
    set fy [interpolateForce grid $r 1]
    set fz [interpolateForce grid $r 2]
    foreach {x y z} $r {break}

    puts $out "$x $y $z $fx $fy $fz"
}
close $out
puts "Done."



