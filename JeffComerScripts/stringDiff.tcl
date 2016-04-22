# Diff two strings.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 2} {
    puts "Usage: tclsh stringDiff.tcl string0 string1"
    exit
}

set s0 [lindex $argv 0]
set s1 [lindex $argv 1]

proc mask {s0 s1} {
    set n0 [string length $s0]
    set n1 [string length $s1]
    set ret ""

    if {$n0 >= $n1} {
	set n [expr $n0]
    } else {
	set n [expr $n1]
    }

    for {set i 0} {$i < $n} {incr i} {
	set c0 [string index $s0 $i]
	set c1 [string index $s1 $i]

	if {[string equal $c0 $c1]} {
	    set ret "${ret}."
	} else {
	    set ret "${ret}x"
	}
    }
    return $ret
}

proc diff {s0 s1} {
    set n0 [string length $s0]
    set n1 [string length $s1]
    
    if {$n0 >= $n1} {
	set n [expr $n0]
    } else {
	set n [expr $n1]
    }

    for {set i 0} {$i < $n} {incr i} {
	set c0 [string index $s0 $i]
	set c1 [string index $s1 $i]

	if {![string equal $c0 $c1]} {
	    puts "c${i} $c0 $c1"
	}
    }
}

puts "Length of string0: [string length $s0]"
puts "Length of string1: [string length $s1]"
puts "Mask: \n[mask $s0 $s1]"
puts "\nDiff: "
diff $s0 $s1

exit



