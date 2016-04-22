# Compute the rmsd of conformations.
# Use with: tclsh rmsdConformations.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set helixLength 3.4
set name conform_ss24_str2
set dispPeriod 100
#Input:
set referenceFile ref.dat
set inFile $name.dat
#Output:
set outFile ref_rmsd_${name}.dat

source vector.tcl

proc rmsd {pos0 pos} {
    set n [llength $pos0]
    set devSum 0.0
    foreach r0 $pos0 r $pos {
	set d [vecLength2 [vecSub $r $r0]]
	set devSum [expr $devSum + $d]
    }
    return [expr sqrt($devSum/$n)]
}

# Read the reference.
set pos0 {}
set in [open $referenceFile r]
foreach line [split [read $in] "\n"] {
    if {[string length $line] < 1} {continue}
    set tok [split [string trim $line]]
    
    if {![string match "!*" $line] && [llength $tok] >= 3} {
	lappend pos0 [lrange $tok 0 2]
    }
}
set n0 [llength $pos0]
close $in

# Read the points in the file.
puts "rmsd"
set nSteps 0
set out [open $outFile w]
set in [open $inFile r]
set pos {}
foreach line [split [read $in] "\n"] {
    if {[string length $line] < 1} {continue}
    set tok [split [string trim $line]]

    if {[string match "!*" $line]} {
	set n [llength $pos]
	if {$n == $n0} {
	    set d [rmsd $pos0 $pos]
	    puts $out "$nSteps $d"
	    if {$nSteps%$dispPeriod == 0} {puts "$nSteps $d"}
	} elseif {$n != 0} {
	    puts "Warning: Reference has $n0 nodes while we read $n nodes."
	}
	
	set pos {}
	incr nSteps
    } elseif {[llength $tok] >= 3} {
	lappend pos [lrange $tok 0 2]
    }
}

# Compute the last one.
if {$n == $n0} {
    set d [rmsd $pos0 $pos]
    puts $out "$nSteps $d"
    if {$nSteps%20 == 0} {puts "$nSteps $d"}
} elseif {$n != 0} {
    puts "Warning: Reference has $n0 nodes while we read $n nodes."
}

puts $out "$nSteps $d"
puts "Computed the rmsd for $nSteps steps."
close $in
close $out
exit



