# Make a sharp grid based on the chosen atoms.
# Use with: tclsh sharpPhantomGrid.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set helixLength 3.4
set name conform_ss0.8_ds2.2_run0
#Input:
set inFile $name.dat
#Output:
set outFile loop_${name}.dat

source vector.tcl
set out [open $outFile w]
set r0 [list 0 0 [expr -$helixLength]]

# Read the points in the file.
puts "Loop-first?"
set nSteps 0
set in [open $inFile r]
set loopFirst 1
set basis [matRandomRot]
set rp0 [vecTransform $basis $r0]
foreach line [split [read $in] "\n"] {
    if {[string length $line] < 1} {continue}
    set tok [split [string trim $line]]

    if {[string match "!*" $line]} {
	set basis [matRandomRot]
	set rp0 [vecTransform $basis $r0]
	puts $out "[expr $nSteps-1] $loopFirst"
	#puts "[expr $nSteps-1] $loopFirst"
	incr nSteps
	set loopFirst 1
    } elseif {[llength $tok] >= 3} {
	set r [list [lindex $tok 0] [lindex $tok 1] [lindex $tok 2]]
	set rp [vecTransform $basis $r]
	set z [lindex $rp 2]
	set z0 [lindex $rp0 2]
	if {$z < 0.0 && $z < $z0} {
	    set loopFirst 0
	}
    }
}
puts "Done."
close $in
close $out
exit



