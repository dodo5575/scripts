# Complement a DNA sequence.
# Runs 5' to 3', so the result is reversed.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 1} {
    puts "Usage: tclsh complement.tcl sequence"
    exit
}

# Complement a sequence.
proc complement {s} {
    set n [string length $s]
    set ret ""

    for {set i 0} {$i < $n} {incr i} {
	set c [string index $s $i]
	
	if {[string equal $c A]} {
	    set c1 T
	} elseif {[string equal $c T]} {
	    set c1 A
	} elseif {[string equal $c G]} {
	    set c1 C
	} elseif {[string equal $c C]} {
	    set c1 G
	} else {
	    puts "Warning: Unrecognized nucleotide $c!"
	    set c1 $c
	}
	
	set ret ${c1}${ret}
    }
    return $ret
}

puts [complement [lindex $argv 0]]



