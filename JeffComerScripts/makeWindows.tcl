# mapWham.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set startInd 300
set posXList {-10 -5 0 5}
set posYList {-10 -5 0 5}
set posZList {13 14 15 16 17 18 19 20}
set spring {0.1 0.1 2.0}
# Output:
set outFile inter_pos.txt

set out [open $outFile w]

set j 0
foreach z $posZList {
    foreach y $posYList {
	foreach x $posXList {
	    puts $out "[expr $startInd+$j] $x $y $z $spring"
	    incr j
	}
    }
}

close $out
