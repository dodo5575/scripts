# mapWham.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set startInd 0
set posXList {-5 -2.5 0 2.5 5}
set posYList {-5 -2.5 0 2.5 5}
set posZList {-10 -5 0 5 10}
set spring {4.0 4.0 4.0}
set ionRadius 2.2
set excludeFile dna_at.txt
# Output:88
set outFile diff_window.txt

source $env(HOME)/scripts/vector.tcl

# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}

	set tok [concat $line]
	lappend r $tok
    }

    close $in
    return $r
}

proc checkPos {posRadList refPos refRad} {
    foreach pos $posRadList {
	set p [lrange $pos 0 2]
	set rad [lindex $pos 3]

	set d [vecLength [vecSub $p $refPos]]
	if {$d < [expr {$rad + $refRad}]} {
	    return 0
	}
    }
    return 1
}

# Read exclusion points.
set exclude [readData $excludeFile]
set out [open $outFile w]

set j 0
foreach z $posZList {
    foreach y $posYList {
	foreach x $posXList {
	    if {[checkPos $exclude [list $x $y $z] $ionRadius]} {
		puts $out "[expr {$startInd+$j}] $x $y $z $spring"
		incr j
	    }
	}
    }
}

puts "Made $j windows."
close $out
